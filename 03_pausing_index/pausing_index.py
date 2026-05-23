#!/usr/bin/env python3
"""
Pausing Index analysis for Chd1 depletion PRO-seq data.
Calculates PI = promoter density / body density for each gene,
compares across 4 conditions, checks replicate quality.
"""

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

PROMOTER_WINDOW = 150  # TSS to +150bp
MIN_GENE_LENGTH = 1000
MIN_READS = 5

SAMPLES = {
    'OsTir1_Ctrl_1':  'dedup/SRR15904219_dedup.bam',
    'OsTir1_Ctrl_2':  'dedup/SRR15904221_dedup.bam',
    'OsTir1_Auxin_1': 'dedup/SRR15904220_dedup.bam',
    'OsTir1_Auxin_2': 'dedup/SRR15904222_dedup.bam',
    'Chd1_Ctrl_1':    'dedup/SRR15904237_dedup.bam',
    'Chd1_Ctrl_2':    'dedup/SRR15904239_dedup.bam',
    'Chd1_Auxin_1':   'dedup/SRR15904238_dedup.bam',
    'Chd1_Auxin_2':   'dedup/SRR15904240_dedup.bam',
}

BED_PLUS = 'annotation/genes_dm6_plus.bed'
BED_MINUS = 'annotation/genes_dm6_minus.bed'
OUT = 'pausing_index'
os.makedirs(OUT, exist_ok=True)


def load_genes(bed_file, strand):
    """Read BED file, keep genes >= 1kb on main chromosomes."""
    genes = []
    valid_chroms = {'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX'}
    with open(bed_file) as f:
        for line in f:
            chrom, start, end, name = line.strip().split('\t')[:4]
            start, end = int(start), int(end)
            if end - start >= MIN_GENE_LENGTH and chrom in valid_chroms:
                genes.append({'chrom': chrom, 'start': start, 'end': end,
                              'name': name, 'strand': strand})
    return genes


def calc_pi(bam_file, genes):
    """
    Calculate Pausing Index for each gene.
    PI = (reads in promoter / promoter length in kb) / (reads in body / body length in kb)
    PRO-seq is reverse-stranded: for + genes, reverse reads are sense.
    """
    bam = pysam.AlignmentFile(bam_file, 'rb')
    results = []
    for g in genes:
        chrom, start, end = g['chrom'], g['start'], g['end']
        strand = g['strand']

        if strand == '+':
            tss = start
            prom = sum(1 for r in bam.fetch(chrom, tss, tss + PROMOTER_WINDOW) if r.is_reverse)
            body = sum(1 for r in bam.fetch(chrom, tss + PROMOTER_WINDOW, end) if r.is_reverse)
        else:
            tss = end
            prom = sum(1 for r in bam.fetch(chrom, tss - PROMOTER_WINDOW, tss) if not r.is_reverse)
            body = sum(1 for r in bam.fetch(chrom, start, tss - PROMOTER_WINDOW) if not r.is_reverse)

        body_len = abs(end - start - PROMOTER_WINDOW)
        if body_len > 0 and prom >= MIN_READS and body >= MIN_READS:
            prom_density = prom / (PROMOTER_WINDOW / 1000)
            body_density = body / (body_len / 1000)
            results.append({'gene': g['name'], 'chrom': chrom,
                            'strand': strand, 'PI': prom_density / body_density})
    bam.close()
    return pd.DataFrame(results)


def merge_replicates(df1, df2):
    """Average PI between two replicates."""
    m = df1[['gene', 'chrom', 'strand', 'PI']].merge(
        df2[['gene', 'PI']], on='gene', suffixes=('_r1', '_r2'))
    m['PI'] = (m['PI_r1'] + m['PI_r2']) / 2
    return m


# load genes
print("Loading genes...")
all_genes = load_genes(BED_PLUS, '+') + load_genes(BED_MINUS, '-')
print(f"  {len(all_genes)} genes")

# calculate PI per sample
print("\nCalculating PI...")
sample_pi = {}
for name, bam_path in SAMPLES.items():
    print(f"  {name}...")
    sample_pi[name] = calc_pi(bam_path, all_genes)
    print(f"    n={len(sample_pi[name])}, median PI={sample_pi[name]['PI'].median():.2f}")

# merge replicates
conditions = {}
for label, r1, r2 in [
    ('OsTir1_Ctrl',  'OsTir1_Ctrl_1',  'OsTir1_Ctrl_2'),
    ('OsTir1_Auxin', 'OsTir1_Auxin_1', 'OsTir1_Auxin_2'),
    ('Chd1_Ctrl',    'Chd1_Ctrl_1',    'Chd1_Ctrl_2'),
    ('Chd1_Auxin',   'Chd1_Auxin_1',   'Chd1_Auxin_2'),
]:
    conditions[label] = merge_replicates(sample_pi[r1], sample_pi[r2])

# keep only genes in all 4 conditions
common_genes = set(conditions['OsTir1_Ctrl']['gene'])
for df in conditions.values():
    common_genes &= set(df['gene'])
for key in conditions:
    conditions[key] = conditions[key][conditions[key]['gene'].isin(common_genes)]

print(f"\nCommon genes: {len(common_genes)}")
for key, df in conditions.items():
    print(f"  {key}: median PI = {df['PI'].median():.2f}")


# pairwise comparisons
print("\nPairwise tests (Wilcoxon signed-rank):")
for c1, c2, desc in [
    ('OsTir1_Ctrl',  'OsTir1_Auxin', 'auxin effect'),
    ('OsTir1_Auxin', 'Chd1_Auxin',   'Chd1 depletion'),
    ('Chd1_Ctrl',    'Chd1_Auxin',   'auxin on Chd1 line'),
]:
    m = conditions[c1][['gene', 'PI']].merge(
        conditions[c2][['gene', 'PI']], on='gene', suffixes=('_1', '_2'))
    _, pval = stats.wilcoxon(m['PI_1'], m['PI_2'])
    fc = m['PI_2'].median() / m['PI_1'].median()
    print(f"  {desc}: FC={fc:.2f}x, p={pval:.2e}")


# --- Plot 1: replicate correlation ---
print("\nPlotting replicate correlations...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, (r1, r2, title) in zip(axes.flat, [
    ('OsTir1_Ctrl_1',  'OsTir1_Ctrl_2',  'OsTir1_Ctrl'),
    ('OsTir1_Auxin_1', 'OsTir1_Auxin_2', 'OsTir1_Auxin'),
    ('Chd1_Ctrl_1',    'Chd1_Ctrl_2',    'Chd1_Ctrl'),
    ('Chd1_Auxin_1',   'Chd1_Auxin_2',   'Chd1_Auxin'),
]):
    m = sample_pi[r1][['gene', 'PI']].merge(
        sample_pi[r2][['gene', 'PI']], on='gene', suffixes=('_1', '_2'))
    x = np.log2(m['PI_1'])
    y = np.log2(m['PI_2'])
    r, _ = stats.pearsonr(x, y)
    ax.scatter(x, y, alpha=0.3, s=3, color='steelblue')
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'r--', alpha=0.5)
    ax.set_title(f'{title} (r={r:.3f})')
    ax.set_xlabel('log2(PI) rep1')
    ax.set_ylabel('log2(PI) rep2')
plt.tight_layout()
plt.savefig(f'{OUT}/PI_replicate_correlation.png', dpi=150)
plt.close()


# --- Plot 2: fold change histograms ---
print("Plotting fold change histograms...")
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, (c1, c2, title) in zip(axes, [
    ('OsTir1_Ctrl',  'OsTir1_Auxin', 'Auxin effect'),
    ('OsTir1_Ctrl',  'Chd1_Ctrl',    'Chd1 ctrl vs OsTir1 ctrl'),
    ('OsTir1_Auxin', 'Chd1_Auxin',   'Chd1 depletion'),
]):
    m = conditions[c1][['gene', 'PI']].merge(
        conditions[c2][['gene', 'PI']], on='gene', suffixes=('_1', '_2'))
    log2fc = np.log2(m['PI_2'] / m['PI_1'])
    ax.hist(log2fc.clip(-5, 5), bins=60, color='steelblue', alpha=0.7, edgecolor='white')
    ax.axvline(0, color='red', ls='--', lw=1.5)
    ax.axvline(log2fc.median(), color='orange', lw=1.5, label=f'Median={log2fc.median():.3f}')
    n_up = (log2fc > 0).sum()
    n_down = (log2fc < 0).sum()
    ax.text(0.05, 0.95, f'Up: {n_up}\nDown: {n_down}', transform=ax.transAxes,
            fontsize=9, va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.set_xlabel(f'log2({c2} / {c1})')
    ax.set_ylabel('Number of genes')
    ax.set_title(title)
    ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(f'{OUT}/PI_foldchange_histograms.png', dpi=150)
plt.close()


# --- Plot 3: X chromosome vs autosomes ---
print("Comparing X vs autosomes...")
merged_xa = conditions['OsTir1_Auxin'][['gene', 'chrom', 'PI']].merge(
    conditions['Chd1_Auxin'][['gene', 'PI']], on='gene', suffixes=('_ctrl', '_mut'))
merged_xa['log2FC'] = np.log2(merged_xa['PI_mut'] / merged_xa['PI_ctrl'])

chrX = merged_xa[merged_xa['chrom'] == 'chrX']
auto = merged_xa[merged_xa['chrom'] != 'chrX']
_, p_xa = stats.mannwhitneyu(chrX['log2FC'], auto['log2FC'])
print(f"  X ({len(chrX)} genes): median log2FC = {chrX['log2FC'].median():.3f}")
print(f"  Autosomes ({len(auto)} genes): median log2FC = {auto['log2FC'].median():.3f}")
print(f"  Mann-Whitney p = {p_xa:.2e}")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (data, title) in zip(axes, [(chrX, 'X chromosome'), (auto, 'Autosomes')]):
    bp = ax.boxplot([np.log2(data['PI_ctrl']), np.log2(data['PI_mut'])],
                    tick_labels=['OsTir1 Auxin', 'Chd1 Auxin'],
                    patch_artist=True, showfliers=False, widths=0.5)
    bp['boxes'][0].set_facecolor('#2166AC'); bp['boxes'][0].set_alpha(0.7)
    bp['boxes'][1].set_facecolor('#D6604D'); bp['boxes'][1].set_alpha(0.7)
    _, p = stats.wilcoxon(data['PI_ctrl'], data['PI_mut'])
    ax.set_ylabel('log2(Pausing Index)')
    ax.set_title(f'{title} (n={len(data)}, p={p:.2e})')
plt.tight_layout()
plt.savefig(f'{OUT}/PI_chrX_vs_autosomes.png', dpi=150)
plt.close()


# --- Plot 4: presentation figure (violin + box) ---
print("Generating presentation figure...")

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.titleweight': 'bold',
    'axes.linewidth': 0.8, 'axes.spines.top': False, 'axes.spines.right': False,
})

colors = {'OsTir1_Ctrl': '#5E81AC', 'OsTir1_Auxin': '#81A1C1',
          'Chd1_Ctrl': '#A3BE8C', 'Chd1_Auxin': '#BF616A'}
cond_order = ['OsTir1_Ctrl', 'OsTir1_Auxin', 'Chd1_Ctrl', 'Chd1_Auxin']
data_pi = [np.log2(conditions[c]['PI'].values) for c in cond_order]

# stats for brackets
m_main = conditions['OsTir1_Auxin'][['gene', 'PI']].merge(
    conditions['Chd1_Auxin'][['gene', 'PI']], on='gene', suffixes=('_1', '_2'))
_, p_main = stats.wilcoxon(m_main['PI_1'], m_main['PI_2'])
fc_main = conditions['Chd1_Auxin']['PI'].median() / conditions['OsTir1_Auxin']['PI'].median()

fig, ax = plt.subplots(figsize=(5, 4.5))

# violin
vp = ax.violinplot(data_pi, positions=[1, 2, 3, 4], showextrema=False, widths=0.7)
for i, body in enumerate(vp['bodies']):
    body.set_facecolor(list(colors.values())[i])
    body.set_alpha(0.3)
    body.set_edgecolor('none')

# box on top
bp = ax.boxplot(data_pi, positions=[1, 2, 3, 4], patch_artist=True,
                showfliers=False, widths=0.25,
                medianprops=dict(color='white', linewidth=1.8),
                whiskerprops=dict(linewidth=0.8, color='#555555'),
                capprops=dict(linewidth=0.8, color='#555555'))
for patch, col in zip(bp['boxes'], colors.values()):
    patch.set_facecolor(col)
    patch.set_alpha(0.85)
    patch.set_edgecolor('#333333')
    patch.set_linewidth(0.7)

# median labels (actual PI values, not log2)
for i, c in enumerate(cond_order):
    med_log = np.median(data_pi[i])
    med_real = conditions[c]['PI'].median()
    ax.text(i + 1, med_log + 0.15, f'{med_real:.1f}', ha='center', va='bottom',
            fontsize=7, color='#555555', style='italic')

# significance brackets
y_top = max(np.percentile(d, 75) + 1.5 * (np.percentile(d, 75) - np.percentile(d, 25))
            for d in data_pi)

# main comparison
ax.plot([2, 2, 4, 4], [y_top + 0.4, y_top + 0.52, y_top + 0.52, y_top + 0.4],
        lw=0.8, c='#333333')
ax.text(3, y_top + 0.57, f'***  FC={fc_main:.1f}x',
        ha='center', va='bottom', fontsize=7.5, color='#333333')

# control comparison
ax.plot([1, 1, 2, 2], [y_top - 0.1, y_top, y_top, y_top - 0.1],
        lw=0.8, c='#333333')
ax.text(1.5, y_top + 0.03, 'n.s.', ha='center', va='bottom',
        fontsize=7.5, color='#888888')

# background shading
ax.axvspan(0.5, 2.5, alpha=0.08, color='#5E81AC', zorder=0)
ax.axvspan(2.5, 4.5, alpha=0.08, color='#BF616A', zorder=0)
ax.text(1.5, ax.get_ylim()[0] + 0.15, 'Controls', ha='center',
        fontsize=7, color='#5E81AC', fontweight='bold', alpha=0.7)
ax.text(3.5, ax.get_ylim()[0] + 0.15, 'Chd1-AID', ha='center',
        fontsize=7, color='#BF616A', fontweight='bold', alpha=0.7)

ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['OsTir1\nCtrl', 'OsTir1\n+ Auxin', 'Chd1-AID\nCtrl', 'Chd1-AID\n+ Auxin'],
                     fontsize=7.5)
ax.set_ylabel(r'$\log_2$(Pausing Index)')
ax.set_title('Promoter-proximal pausing (PI)')

plt.tight_layout()
plt.savefig(f'{OUT}/PI_figure.png', dpi=300, facecolor='white')
plt.close()


# summary
print("\nSummary:")
for key in cond_order:
    print(f"  {key}: median PI = {conditions[key]['PI'].median():.2f}")
print(f"  FC (Chd1_Auxin / OsTir1_Auxin) = {fc_main:.2f}x, p = {p_main:.2e}")
print(f"  Genes analyzed: {len(common_genes)}")
print(f"\nPlots saved to {OUT}/")
