#!/usr/bin/env python3
"""
Aberrant (antisense) transcription analysis across 4 conditions
Counts sense/antisense reads in gene bodies, compares between conditions
Triple control: OsTir1_Ctrl, OsTir1_Auxin, Chd1_Ctrl, Chd1_Auxin
"""

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os

MIN_GENE_LENGTH = 1000
MIN_SENSE_READS = 10
BODY_OFFSET = 500  # skip promoter region

SAMPLES = {
    "OsTir1_Ctrl_1": "dedup/SRR15904219_dedup.bam",
    "OsTir1_Ctrl_2": "dedup/SRR15904221_dedup.bam",
    "OsTir1_Auxin_1": "dedup/SRR15904220_dedup.bam",
    "OsTir1_Auxin_2": "dedup/SRR15904222_dedup.bam",
    "Chd1_Ctrl_1": "dedup/SRR15904237_dedup.bam",
    "Chd1_Ctrl_2": "dedup/SRR15904239_dedup.bam",
    "Chd1_Auxin_1": "dedup/SRR15904238_dedup.bam",
    "Chd1_Auxin_2": "dedup/SRR15904240_dedup.bam",
}

BED_PLUS = "annotation/genes_dm6_plus.bed"
BED_MINUS = "annotation/genes_dm6_minus.bed"
OUT = "aberrant_transcription"
os.makedirs(OUT, exist_ok=True)


def load_genes(bed_file, strand):
    genes = []
    valid_chroms = {"chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"}
    with open(bed_file) as f:
        for line in f:
            chrom, start, end, name = line.strip().split("\t")[:4]
            start, end = int(start), int(end)
            if end - start >= MIN_GENE_LENGTH and chrom in valid_chroms:
                genes.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "name": name,
                        "strand": strand,
                    }
                )
    return genes


def count_sense_antisense(bam_file, genes):
    # PRO-seq is reverse-stranded: for + genes, reverse reads = sense
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = []
    for g in genes:
        chrom, strand = g["chrom"], g["strand"]
        if strand == "+":
            body_start = g["start"] + BODY_OFFSET
            body_end = g["end"]
        else:
            body_start = g["start"]
            body_end = g["end"] - BODY_OFFSET
        if body_end <= body_start:
            continue

        sense, antisense = 0, 0
        for read in bam.fetch(chrom, body_start, body_end):
            if strand == "+":
                if read.is_reverse:
                    sense += 1
                else:
                    antisense += 1
            else:
                if not read.is_reverse:
                    sense += 1
                else:
                    antisense += 1

        if sense >= MIN_SENSE_READS:
            results.append(
                {
                    "gene": g["name"],
                    "chrom": chrom,
                    "strand": strand,
                    "sense": sense,
                    "antisense": antisense,
                    "antisense_ratio": antisense / sense,
                }
            )
    bam.close()
    return pd.DataFrame(results)


def merge_replicates(df1, df2, label):
    merged = df1[["gene", "chrom", "strand", "antisense_ratio"]].merge(
        df2[["gene", "antisense_ratio"]], on="gene", suffixes=("_r1", "_r2")
    )
    merged["antisense_ratio"] = (
        merged["antisense_ratio_r1"] + merged["antisense_ratio_r2"]
    ) / 2
    merged["condition"] = label
    return merged[["gene", "chrom", "strand", "antisense_ratio", "condition"]]


# load genes
print("Loading genes...")
all_genes = load_genes(BED_PLUS, "+") + load_genes(BED_MINUS, "-")
print(f"  {len(all_genes)} genes loaded")

# count reads for each sample
print("\nCounting reads...")
sample_data = {}
for name, bam_path in SAMPLES.items():
    print(f"  {name}...")
    sample_data[name] = count_sense_antisense(bam_path, all_genes)
    med = sample_data[name]["antisense_ratio"].median()
    print(f"    n={len(sample_data[name])}, median AS ratio={med:.4f}")

# merge replicates
print("\nMerging replicates...")
conditions = {}
for label, r1, r2 in [
    ("OsTir1_Ctrl", "OsTir1_Ctrl_1", "OsTir1_Ctrl_2"),
    ("OsTir1_Auxin", "OsTir1_Auxin_1", "OsTir1_Auxin_2"),
    ("Chd1_Ctrl", "Chd1_Ctrl_1", "Chd1_Ctrl_2"),
    ("Chd1_Auxin", "Chd1_Auxin_1", "Chd1_Auxin_2"),
]:
    conditions[label] = merge_replicates(sample_data[r1], sample_data[r2], label)

# only keep genes present in all 4 conditions
common_genes = set(conditions["OsTir1_Ctrl"]["gene"])
for df in conditions.values():
    common_genes &= set(df["gene"])
for key in conditions:
    conditions[key] = conditions[key][conditions[key]["gene"].isin(common_genes)]

print(f"  Common genes: {len(common_genes)}")
for key, df in conditions.items():
    print(f"  {key}: median = {df['antisense_ratio'].median():.4f}")


# pairwise comparisons
print("\nPairwise comparisons (Wilcoxon signed-rank):")
for c1, c2, desc in [
    ("OsTir1_Ctrl", "OsTir1_Auxin", "auxin effect"),
    ("OsTir1_Auxin", "Chd1_Auxin", "Chd1 depletion"),
    ("Chd1_Ctrl", "Chd1_Auxin", "auxin on Chd1 line"),
]:
    m = conditions[c1][["gene", "antisense_ratio"]].merge(
        conditions[c2][["gene", "antisense_ratio"]], on="gene", suffixes=("_1", "_2")
    )
    _, pval = stats.wilcoxon(m["antisense_ratio_1"], m["antisense_ratio_2"])
    fc = m["antisense_ratio_2"].median() / m["antisense_ratio_1"].median()
    n_up = (m["antisense_ratio_2"] > m["antisense_ratio_1"]).sum()
    print(f"  {desc}: FC={fc:.2f}x, genes up={n_up}/{len(m)}, p={pval:.2e}")


# replicate correlation
print("\nPlotting replicate correlations...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, (r1, r2, title) in zip(
    axes.flat,
    [
        ("OsTir1_Ctrl_1", "OsTir1_Ctrl_2", "OsTir1_Ctrl"),
        ("OsTir1_Auxin_1", "OsTir1_Auxin_2", "OsTir1_Auxin"),
        ("Chd1_Ctrl_1", "Chd1_Ctrl_2", "Chd1_Ctrl"),
        ("Chd1_Auxin_1", "Chd1_Auxin_2", "Chd1_Auxin"),
    ],
):
    m = sample_data[r1][["gene", "antisense_ratio"]].merge(
        sample_data[r2][["gene", "antisense_ratio"]], on="gene", suffixes=("_1", "_2")
    )
    x = np.log2(m["antisense_ratio_1"] + 0.001)
    y = np.log2(m["antisense_ratio_2"] + 0.001)
    r, _ = stats.pearsonr(x, y)
    ax.scatter(x, y, alpha=0.3, s=3, color="steelblue")
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, "r--", alpha=0.5)
    ax.set_title(f"{title} (r={r:.3f})")
    ax.set_xlabel("log2(AS ratio) rep1")
    ax.set_ylabel("log2(AS ratio) rep2")

plt.tight_layout()
plt.savefig(f"{OUT}/replicate_correlation.png", dpi=150)
plt.close()


# boxplot
cond_order = ["OsTir1_Ctrl", "OsTir1_Auxin", "Chd1_Ctrl", "Chd1_Auxin"]
plot_data = [np.log2(conditions[c]["antisense_ratio"] + 0.001) for c in cond_order]

fig, ax = plt.subplots(figsize=(8, 6))
bp = ax.boxplot(
    plot_data,
    patch_artist=True,
    showfliers=False,
    widths=0.6,
    tick_labels=["OsTir1 Ctrl", "OsTir1 Auxin", "Chd1 Ctrl", "Chd1 Auxin"],
)
for patch, col in zip(bp["boxes"], ["#2166AC", "#4393C3", "#1A9641", "#D6604D"]):
    patch.set_facecolor(col)
    patch.set_alpha(0.7)
ax.set_ylabel("log2(Antisense / Sense ratio)")
ax.set_title("Aberrant transcription across conditions")
plt.tight_layout()
plt.savefig(f"{OUT}/AT_boxplot.png", dpi=150)
plt.close()


# fold change histograms
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, (c1, c2, title) in zip(
    axes,
    [
        ("OsTir1_Ctrl", "OsTir1_Auxin", "Auxin effect"),
        ("OsTir1_Ctrl", "Chd1_Ctrl", "Chd1 ctrl vs OsTir1 ctrl"),
        ("OsTir1_Auxin", "Chd1_Auxin", "Chd1 depletion"),
    ],
):
    m = conditions[c1][["gene", "antisense_ratio"]].merge(
        conditions[c2][["gene", "antisense_ratio"]], on="gene", suffixes=("_1", "_2")
    )
    log2fc = np.log2(
        (m["antisense_ratio_2"] + 0.001) / (m["antisense_ratio_1"] + 0.001)
    )
    ax.hist(
        log2fc.clip(-5, 5), bins=60, color="steelblue", alpha=0.7, edgecolor="white"
    )
    ax.axvline(0, color="red", ls="--", lw=1.5)
    ax.axvline(
        log2fc.median(), color="orange", lw=1.5, label=f"Median={log2fc.median():.3f}"
    )
    n_up = (log2fc > 0).sum()
    n_down = (log2fc < 0).sum()
    ax.text(
        0.05,
        0.95,
        f"Up: {n_up}\nDown: {n_down}",
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )
    ax.set_xlabel(f"log2({c2} / {c1})")
    ax.set_ylabel("Number of genes")
    ax.set_title(title)
    ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(f"{OUT}/AT_foldchange_histograms.png", dpi=150)
plt.close()


# scatter: main comparison
m = conditions["OsTir1_Auxin"][["gene", "antisense_ratio"]].merge(
    conditions["Chd1_Auxin"][["gene", "antisense_ratio"]],
    on="gene",
    suffixes=("_ctrl", "_mut"),
)
fig, ax = plt.subplots(figsize=(7, 7))
ax.scatter(
    np.log2(m["antisense_ratio_ctrl"] + 0.001),
    np.log2(m["antisense_ratio_mut"] + 0.001),
    alpha=0.3,
    s=5,
    color="steelblue",
)
lims = [
    min(ax.get_xlim()[0], ax.get_ylim()[0]),
    max(ax.get_xlim()[1], ax.get_ylim()[1]),
]
ax.plot(lims, lims, "r--", alpha=0.5, label="x = y")
ax.set_xlabel("log2(AS ratio) OsTir1 Auxin")
ax.set_ylabel("log2(AS ratio) Chd1 Auxin")
ax.set_title("Antisense transcription: gene-by-gene")
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/AT_scatter.png", dpi=150)
plt.close()


# X chromosome vs autosomes
merged = conditions["OsTir1_Auxin"][["gene", "chrom", "antisense_ratio"]].merge(
    conditions["Chd1_Auxin"][["gene", "antisense_ratio"]],
    on="gene",
    suffixes=("_ctrl", "_mut"),
)
merged["log2FC"] = np.log2(
    (merged["antisense_ratio_mut"] + 0.001) / (merged["antisense_ratio_ctrl"] + 0.001)
)

chrX = merged[merged["chrom"] == "chrX"]
auto = merged[merged["chrom"] != "chrX"]
print(
    f"\nX chromosome ({len(chrX)} genes): median log2FC = {chrX['log2FC'].median():.3f}"
)
print(f"Autosomes ({len(auto)} genes): median log2FC = {auto['log2FC'].median():.3f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (data, title) in zip(axes, [(chrX, "X chromosome"), (auto, "Autosomes")]):
    bp = ax.boxplot(
        [
            np.log2(data["antisense_ratio_ctrl"] + 0.001),
            np.log2(data["antisense_ratio_mut"] + 0.001),
        ],
        tick_labels=["OsTir1 Auxin", "Chd1 Auxin"],
        patch_artist=True,
        showfliers=False,
        widths=0.5,
    )
    bp["boxes"][0].set_facecolor("#2166AC")
    bp["boxes"][0].set_alpha(0.7)
    bp["boxes"][1].set_facecolor("#D6604D")
    bp["boxes"][1].set_alpha(0.7)
    _, p = stats.wilcoxon(data["antisense_ratio_ctrl"], data["antisense_ratio_mut"])
    ax.set_ylabel("log2(Antisense / Sense)")
    ax.set_title(f"{title} (n={len(data)}, p={p:.2e})")
plt.tight_layout()
plt.savefig(f"{OUT}/AT_chrX_vs_autosomes.png", dpi=150)
plt.close()


# top genes
print("\nTop 20 genes with highest antisense increase:")
top20 = merged.nlargest(20, "log2FC")
print(f"{'Gene':<15} {'Chrom':<8} {'AS_ctrl':>10} {'AS_mut':>10} {'log2FC':>8}")
for _, row in top20.iterrows():
    print(
        f"{row['gene']:<15} {row['chrom']:<8} "
        f"{row['antisense_ratio_ctrl']:>10.4f} "
        f"{row['antisense_ratio_mut']:>10.4f} "
        f"{row['log2FC']:>8.2f}"
    )

# summary
print("\nMedian AS ratio per condition:")
for key in cond_order:
    med = conditions[key]["antisense_ratio"].median()
    print(f"  {key}: {med:.4f} ({med*100:.1f}%)")
print(f"\nPlots saved to {OUT}/")

# presentation figure (violin + box, styled)
print("Generating main figure...")

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

at_colors = ["#5E81AC", "#81A1C1", "#A3BE8C", "#BF616A"]
at_labels = ["OsTir1\nCtrl", "OsTir1\n+ Auxin", "Chd1-AID\nCtrl", "Chd1-AID\n+ Auxin"]
data_styled = [
    np.log2(conditions[c]["antisense_ratio"].values + 0.001) for c in cond_order
]

fig, ax = plt.subplots(figsize=(5, 4.5))
ax.set_ylim(-8, 9)

vp = ax.violinplot(data_styled, positions=[1, 2, 3, 4], showextrema=False, widths=0.7)
for i, body in enumerate(vp["bodies"]):
    body.set_facecolor(at_colors[i])
    body.set_alpha(0.35)
    body.set_edgecolor("none")

bp = ax.boxplot(
    data_styled,
    positions=[1, 2, 3, 4],
    patch_artist=True,
    showfliers=False,
    widths=0.25,
    medianprops=dict(color="white", linewidth=1.8),
    whiskerprops=dict(linewidth=0.8, color="#555555"),
    capprops=dict(linewidth=0.8, color="#555555"),
)
for patch, col in zip(bp["boxes"], at_colors):
    patch.set_facecolor(col)
    patch.set_alpha(0.85)
    patch.set_edgecolor("#333333")
    patch.set_linewidth(0.7)

for i, d in enumerate(data_styled):
    med = np.median(d)
    ax.text(
        i + 1,
        med - 0.35,
        f"{med:.2f}",
        ha="center",
        va="top",
        fontsize=6.5,
        color="#555555",
        style="italic",
    )

# stats for brackets
m_main = conditions["OsTir1_Auxin"][["gene", "antisense_ratio"]].merge(
    conditions["Chd1_Auxin"][["gene", "antisense_ratio"]],
    on="gene",
    suffixes=("_1", "_2"),
)
_, p_main = stats.wilcoxon(m_main["antisense_ratio_1"], m_main["antisense_ratio_2"])
fc_at = (
    conditions["Chd1_Auxin"]["antisense_ratio"].median()
    / conditions["OsTir1_Auxin"]["antisense_ratio"].median()
)

y_top = max(
    np.percentile(d, 75) + 1.5 * (np.percentile(d, 75) - np.percentile(d, 25))
    for d in data_styled
)

# main bracket
ax.plot(
    [2, 2, 4, 4],
    [y_top + 0.4, y_top + 0.52, y_top + 0.52, y_top + 0.4],
    lw=0.8,
    c="#333333",
)
ax.text(
    3,
    y_top + 0.57,
    f"***  FC={fc_at:.1f}x",
    ha="center",
    va="bottom",
    fontsize=7.5,
    color="#333333",
)

# control bracket
ax.plot([1, 1, 2, 2], [y_top - 0.1, y_top, y_top, y_top - 0.1], lw=0.8, c="#333333")
ax.text(
    1.5, y_top + 0.03, "n.s.", ha="center", va="bottom", fontsize=7.5, color="#888888"
)

# background
ax.axvspan(0.5, 2.5, alpha=0.08, color="#5E81AC", zorder=0)
ax.axvspan(2.5, 4.5, alpha=0.08, color="#BF616A", zorder=0)
ax.text(
    1.5,
    ax.get_ylim()[0] + 0.15,
    "Controls",
    ha="center",
    fontsize=7,
    color="#5E81AC",
    fontweight="bold",
    alpha=0.7,
)
ax.text(
    3.5,
    ax.get_ylim()[0] + 0.15,
    "Chd1-AID",
    ha="center",
    fontsize=7,
    color="#BF616A",
    fontweight="bold",
    alpha=0.7,
)

ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(at_labels, fontsize=7.5)
ax.set_ylabel(r"$\log_2$(Antisense / Sense ratio)")
ax.set_title("Aberrant (antisense) transcription")

plt.tight_layout()
fig.savefig(f"{OUT}/AT_figure.png", dpi=300, facecolor="white")
plt.close()
