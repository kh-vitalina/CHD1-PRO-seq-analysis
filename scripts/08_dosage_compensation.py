#!/usr/bin/env python3
"""
Dosage compensation analysis for Chd1 depletion PRO-seq data.
Compares transcription levels on X chromosome vs autosomes.
In Drosophila males (and S2 cells), MSL complex upregulates X ~2x.
"""

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

MIN_GENE_LENGTH = 1000
MIN_READS = 10

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
OUT = 'dosage_compensation'
os.makedirs(OUT, exist_ok=True)


def load_genes(bed_file, strand):
    """Read BED, keep genes >= 1kb on main chromosomes."""
    genes = []
    valid_chroms = {'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX'}
    with open(bed_file) as f:
        for line in f:
            chrom, start, end, name = line.strip().split('\t')[:4]
            start, end = int(start), int(end)
            if end - start >= MIN_GENE_LENGTH and chrom in valid_chroms:
                genes.append({'chrom': chrom, 'start': start, 'end': end,
                              'name': name, 'strand': strand, 'length': end - start})
    return genes


def count_reads_rpk(bam_file, genes):
    """Count sense-strand reads, normalize by gene length (RPK)."""
    bam = pysam.AlignmentFile(bam_file, 'rb')
    results = []
    for g in genes:
        chrom, start, end = g['chrom'], g['start'], g['end']
        strand, length = g['strand'], g['length']
        count = 0
        for read in bam.fetch(chrom, start, end):
            if strand == '+' and read.is_reverse:
                count += 1
            elif strand == '-' and not read.is_reverse:
                count += 1
        if count >= MIN_READS:
            results.append({
                'gene': g['name'], 'chrom': chrom, 'strand': strand,
                'length': length, 'reads': count, 'rpk': count / (length / 1000)
            })
    bam.close()
    return pd.DataFrame(results)


def merge_replicates(df1, df2):
    """Average RPK between two replicates."""
    m = df1[['gene', 'chrom', 'strand', 'length', 'rpk']].merge(
        df2[['gene', 'rpk']], on='gene', suffixes=('_r1', '_r2'))
    m['rpk'] = (m['rpk_r1'] + m['rpk_r2']) / 2
    return m


# load genes
print("Loading genes...")
all_genes = load_genes(BED_PLUS, '+') + load_genes(BED_MINUS, '-')
n_x = sum(1 for g in all_genes if g['chrom'] == 'chrX')
print(f"  {len(all_genes)} genes ({n_x} on X, {len(all_genes) - n_x} on autosomes)")

# count reads per sample
print("\nCounting reads...")
sample_data = {}
for name, bam_path in SAMPLES.items():
    print(f"  {name}...")
    sample_data[name] = count_reads_rpk(bam_path, all_genes)
    df = sample_data[name]
    x_med = df[df['chrom'] == 'chrX']['rpk'].median()
    a_med = df[df['chrom'] != 'chrX']['rpk'].median()
    print(f"    n={len(df)}, X/A={x_med/a_med:.3f}")

# merge replicates
conditions = {}
for label, r1, r2 in [
    ('OsTir1_Ctrl',  'OsTir1_Ctrl_1',  'OsTir1_Ctrl_2'),
    ('OsTir1_Auxin', 'OsTir1_Auxin_1', 'OsTir1_Auxin_2'),
    ('Chd1_Ctrl',    'Chd1_Ctrl_1',    'Chd1_Ctrl_2'),
    ('Chd1_Auxin',   'Chd1_Auxin_1',   'Chd1_Auxin_2'),
]:
    conditions[label] = merge_replicates(sample_data[r1], sample_data[r2])

# common genes only
common_genes = set(conditions['OsTir1_Ctrl']['gene'])
for df in conditions.values():
    common_genes &= set(df['gene'])
for key in conditions:
    conditions[key] = conditions[key][conditions[key]['gene'].isin(common_genes)]
print(f"\nCommon genes: {len(common_genes)}")


# X/A ratio per condition
print("\nX/A expression ratios:")
xa_ratios = {}
for key in ['OsTir1_Ctrl', 'OsTir1_Auxin', 'Chd1_Ctrl', 'Chd1_Auxin']:
    df = conditions[key]
    x_med = df[df['chrom'] == 'chrX']['rpk'].median()
    a_med = df[df['chrom'] != 'chrX']['rpk'].median()
    xa_ratios[key] = x_med / a_med
    print(f"  {key}: X/A = {xa_ratios[key]:.3f}")


# gene-level fold changes
merged = conditions['OsTir1_Auxin'][['gene', 'chrom', 'rpk']].merge(
    conditions['Chd1_Auxin'][['gene', 'rpk']], on='gene', suffixes=('_ctrl', '_mut'))
merged['log2FC'] = np.log2((merged['rpk_mut'] + 0.1) / (merged['rpk_ctrl'] + 0.1))

chrX = merged[merged['chrom'] == 'chrX']
auto = merged[merged['chrom'] != 'chrX']
_, p_xa = stats.mannwhitneyu(chrX['log2FC'], auto['log2FC'])

print(f"\nFold change upon Chd1 depletion:")
print(f"  X ({len(chrX)} genes): median log2FC = {chrX['log2FC'].median():.3f}")
print(f"  Autosomes ({len(auto)} genes): median log2FC = {auto['log2FC'].median():.3f}")
print(f"  Mann-Whitney p = {p_xa:.2e}")


# per-chromosome breakdown
print("\nPer-chromosome fold changes:")
chrom_order = ['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']
chrom_stats = {}
for chrom in chrom_order:
    sub = merged[merged['chrom'] == chrom]
    _, p = stats.wilcoxon(sub['rpk_ctrl'], sub['rpk_mut']) if len(sub) > 10 else (None, 1.0)
    chrom_stats[chrom] = {'n': len(sub), 'median': sub['log2FC'].median(), 'p': p}
    print(f"  {chrom}: n={len(sub)}, median log2FC={sub['log2FC'].median():.3f}, p={p:.2e}")


# --- Plot 1: X/A ratio barplot ---
print("\nPlotting X/A ratios...")
fig, ax = plt.subplots(figsize=(8, 6))
cond_names = ['OsTir1_Ctrl', 'OsTir1_Auxin', 'Chd1_Ctrl', 'Chd1_Auxin']
ratios = [xa_ratios[c] for c in cond_names]
bar_colors = ['#2166AC', '#4393C3', '#1A9641', '#D6604D']
labels = ['OsTir1 Ctrl', 'OsTir1 Auxin', 'Chd1 Ctrl', 'Chd1 Auxin']

bars = ax.bar(range(len(ratios)), ratios, color=bar_colors, alpha=0.7, edgecolor='black')
ax.set_xticks(range(len(ratios)))
ax.set_xticklabels(labels)
ax.set_ylabel('X / Autosome expression ratio')
ax.set_title('Dosage compensation: X/A ratio')
ax.axhline(1.0, color='grey', ls='--', alpha=0.5)
for bar, val in zip(bars, ratios):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{val:.3f}', ha='center', fontsize=10, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{OUT}/XA_ratio_barplot.png', dpi=150)
plt.close()


# --- Plot 2: fold change boxplot X vs autosomes ---
fig, ax = plt.subplots(figsize=(8, 6))
bp = ax.boxplot([chrX['log2FC'], auto['log2FC']],
                tick_labels=[f'X (n={len(chrX)})', f'Autosomes (n={len(auto)})'],
                patch_artist=True, showfliers=False, widths=0.5)
bp['boxes'][0].set_facecolor('#E41A1C'); bp['boxes'][0].set_alpha(0.7)
bp['boxes'][1].set_facecolor('#377EB8'); bp['boxes'][1].set_alpha(0.7)
ax.axhline(0, color='grey', ls='--', alpha=0.5)
ax.set_ylabel('log2(Chd1 Auxin / OsTir1 Auxin)')
ax.set_title(f'Expression change: X vs Autosomes (p={p_xa:.2e})')
plt.tight_layout()
plt.savefig(f'{OUT}/XA_foldchange_boxplot.png', dpi=150)
plt.close()


# --- Plot 3: per-chromosome barplot ---
fig, ax = plt.subplots(figsize=(8, 5))
medians = [chrom_stats[c]['median'] for c in chrom_order]
chrom_colors = ['#E41A1C'] + ['#377EB8'] * 5
ax.bar(range(len(chrom_order)), medians, color=chrom_colors, alpha=0.7, edgecolor='black')
ax.set_xticks(range(len(chrom_order)))
ax.set_xticklabels(chrom_order)
ax.axhline(0, color='grey', ls='--', alpha=0.5)
ax.set_ylabel('Median log2FC upon Chd1 depletion')
ax.set_title('Transcription change by chromosome')
for i, med in enumerate(medians):
    ax.text(i, med + 0.01, f'{med:.3f}', ha='center', fontsize=9, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{OUT}/per_chromosome_barplot.png', dpi=150)
plt.close()


# --- Plot 4: presentation figure (paired barplot) ---
print("Generating presentation figure...")

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.titleweight': 'bold',
    'axes.linewidth': 0.8, 'axes.spines.top': False, 'axes.spines.right': False,
})

fig, ax = plt.subplots(figsize=(6, 4.5))
chroms = ['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']
chrom_nice = ['X', '2L', '2R', '3L', '3R', '4']
x = np.arange(len(chroms))
width = 0.35

med_ctrl = []
med_dep = []
n_per_chrom = []
for chrom in chroms:
    c = conditions['OsTir1_Auxin']
    d = conditions['Chd1_Auxin']
    med_ctrl.append(c[c['chrom'] == chrom]['rpk'].median())
    med_dep.append(d[d['chrom'] == chrom]['rpk'].median())
    n_per_chrom.append(len(c[c['chrom'] == chrom]))

ax.bar(x - width/2, med_ctrl, width, label='OsTir1 + Auxin',
       color='#5E81AC', alpha=0.8, edgecolor='white', linewidth=1)
ax.bar(x + width/2, med_dep, width, label='Chd1 depleted',
       color='#BF616A', alpha=0.8, edgecolor='white', linewidth=1)

ax.axvspan(-0.5, 0.5, alpha=0.06, color='#C0392B', zorder=0)

for i, n in enumerate(n_per_chrom):
    y = max(med_ctrl[i], med_dep[i])
    ax.text(i, y + 0.5, f'n={n}', ha='center', va='bottom', fontsize=6, color='#777777')

ax.set_xticks(x)
ax.set_xticklabels(chrom_nice, fontsize=9)
ax.set_xlabel('Chromosome')
ax.set_ylabel('Median RPK')
ax.set_title('Transcription level by chromosome')
ax.legend(loc='upper left', fontsize=8, bbox_to_anchor=(0.0, 0.78))

xa_c = med_ctrl[0] / np.median(med_ctrl[1:])
xa_d = med_dep[0] / np.median(med_dep[1:])
ax.text(0.02, 0.97,
        f'X/A ratio:\nControl: {xa_c:.2f}\nChd1 dep: {xa_d:.2f}\np = {p_xa:.1e}',
        transform=ax.transAxes, fontsize=7, ha='left', va='top', color='#C0392B')

plt.tight_layout()
plt.savefig(f'{OUT}/DC_figure.png', dpi=300, facecolor='white')
plt.close()


# summary
print("\nSummary:")
print(f"  X/A ratio: control={xa_ratios['OsTir1_Auxin']:.3f}, "
      f"Chd1 dep={xa_ratios['Chd1_Auxin']:.3f}")
print(f"  X vs Auto p = {p_xa:.2e}")
print(f"  Genes analyzed: {len(common_genes)}")
print(f"\nPlots saved to {OUT}/")
