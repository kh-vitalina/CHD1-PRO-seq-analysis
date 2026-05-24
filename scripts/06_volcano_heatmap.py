#!/usr/bin/env python3
"""
Volcano plot and heatmap for DESeq2 results.
Volcano: Chd1 depletion contrast with gene labels.
Heatmap: top 50 up + 50 down genes, z-score normalized.
Also produces sample correlation heatmap.
"""

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import os

try:
    from adjustText import adjust_text

    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False
    print("adjustText not installed, gene labels may overlap")

TOP_N = 50
MIN_GENE_LENGTH = 1000

SAMPLES_ORDER = [
    ("OsTir1_Ctrl_1", "dedup/SRR15904219_dedup.bam"),
    ("OsTir1_Ctrl_2", "dedup/SRR15904221_dedup.bam"),
    ("OsTir1_Auxin_1", "dedup/SRR15904220_dedup.bam"),
    ("OsTir1_Auxin_2", "dedup/SRR15904222_dedup.bam"),
    ("Chd1_Ctrl_1", "dedup/SRR15904237_dedup.bam"),
    ("Chd1_Ctrl_2", "dedup/SRR15904239_dedup.bam"),
    ("Chd1_Auxin_1", "dedup/SRR15904238_dedup.bam"),
    ("Chd1_Auxin_2", "dedup/SRR15904240_dedup.bam"),
]

BED_PLUS = "annotation/genes_dm6_plus.bed"
BED_MINUS = "annotation/genes_dm6_minus.bed"
OUT = "deseq2_results"
os.makedirs(OUT, exist_ok=True)


def load_genes(bed_file, strand):
    """Read BED, keep genes >= 1kb on main chromosomes."""
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


# load DESeq2 results
deseq2_file = "deseq2_results/DESeq2_Chd1_depletion.csv"
print("Loading DESeq2 results...")
if not os.path.exists(deseq2_file):
    print(f"ERROR: {deseq2_file} not found, run deseq2_analysis.R first")
    exit(1)

deseq2 = pd.read_csv(deseq2_file).dropna(subset=["log2FoldChange", "padj"])

# figure out gene column name
gene_col = "gene"
if gene_col not in deseq2.columns:
    for col in ["gene_name", "name", "Gene", "Name"]:
        if col in deseq2.columns:
            gene_col = col
            break
    else:
        deseq2["gene"] = deseq2.index
        gene_col = "gene"

print(f"  {len(deseq2)} genes, column: '{gene_col}'")


# --- volcano plot ---
print("\nPlotting volcano...")

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

fig, ax = plt.subplots(figsize=(6, 4))

nonsig = deseq2[deseq2["padj"] >= 0.05]
sig = deseq2[(deseq2["padj"] < 0.05) & (deseq2["log2FoldChange"].abs() < 1)]
sig_up = deseq2[(deseq2["padj"] < 0.05) & (deseq2["log2FoldChange"] >= 1)]
sig_down = deseq2[(deseq2["padj"] < 0.05) & (deseq2["log2FoldChange"] <= -1)]

ax.scatter(
    nonsig["log2FoldChange"],
    -np.log10(nonsig["padj"]),
    s=3,
    alpha=0.25,
    color="#D0D0D0",
    edgecolors="none",
    rasterized=True,
)
ax.scatter(
    sig["log2FoldChange"],
    -np.log10(sig["padj"]),
    s=3,
    alpha=0.45,
    color="#F4A582",
    edgecolors="none",
    rasterized=True,
)
ax.scatter(
    sig_up["log2FoldChange"],
    -np.log10(sig_up["padj"]),
    s=6,
    alpha=0.7,
    color="#D73027",
    edgecolors="none",
    rasterized=True,
)
ax.scatter(
    sig_down["log2FoldChange"],
    -np.log10(sig_down["padj"]),
    s=6,
    alpha=0.7,
    color="#4575B4",
    edgecolors="none",
    rasterized=True,
)

ax.axhline(-np.log10(0.05), ls="--", color="#999999", lw=0.5)
ax.axvline(-1, ls="--", color="#999999", lw=0.5)
ax.axvline(1, ls="--", color="#999999", lw=0.5)

ax.text(
    0.96,
    0.96,
    f"↑ {len(sig_up)}",
    transform=ax.transAxes,
    fontsize=9,
    ha="right",
    va="top",
    color="#D73027",
    fontweight="bold",
)
ax.text(
    0.04,
    0.96,
    f"↓ {len(sig_down)}",
    transform=ax.transAxes,
    fontsize=9,
    ha="left",
    va="top",
    color="#4575B4",
    fontweight="bold",
)

# label top genes
top_label = deseq2.copy()
top_label["nlp"] = -np.log10(top_label["padj"].clip(lower=1e-300))
top_label = top_label.nlargest(15, "nlp")

if HAS_ADJUSTTEXT:
    texts = []
    for _, row in top_label.iterrows():
        y_val = min(
            -np.log10(max(row["padj"], 1e-300)),
            ax.get_ylim()[1] - 5 if ax.get_ylim()[1] > 10 else 300,
        )
        texts.append(
            ax.annotate(
                row[gene_col],
                (row["log2FoldChange"], y_val),
                fontsize=6.75,
                color="#333333",
                ha="center",
                va="bottom",
                xytext=(0, 3),
                textcoords="offset points",
            )
        )
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5))
else:
    for _, row in top_label.iterrows():
        y_val = min(-np.log10(max(row["padj"], 1e-300)), 280)
        ax.annotate(
            row[gene_col],
            (row["log2FoldChange"], y_val),
            fontsize=5,
            color="#333333",
            ha="center",
            va="bottom",
            xytext=(0, 3),
            textcoords="offset points",
        )

ax.set_xlabel(r"$\log_2$(Fold Change)")
ax.set_ylabel(r"$-\log_{10}$(adjusted p-value)")
ax.set_title("Differential transcription upon Chd1 depletion")

plt.tight_layout()
fig.savefig(f"{OUT}/volcano_chd1_depletion.png", dpi=300, facecolor="white")
plt.close()

# print top genes
deseq2["nlp"] = -np.log10(deseq2["padj"].clip(lower=1e-300))
top20 = deseq2.nlargest(20, "nlp")
print("\nTop 20 most significant genes:")
print(top20[[gene_col, "log2FoldChange", "padj"]].to_string())


# --- heatmap of top DE genes ---
print("\nPreparing heatmap...")

# select top genes by fold change
sig_genes = deseq2[deseq2["padj"] < 0.05].copy()
sig_genes["abs_lfc"] = sig_genes["log2FoldChange"].abs()
top_up = sig_genes[sig_genes["log2FoldChange"] > 0].nlargest(TOP_N, "abs_lfc")
top_down = sig_genes[sig_genes["log2FoldChange"] < 0].nlargest(TOP_N, "abs_lfc")
top_genes_df = pd.concat([top_up, top_down])
selected = set(top_genes_df[gene_col].tolist())
print(f"  Selected {len(selected)} genes ({len(top_up)} up + {len(top_down)} down)")

# load annotation and match
all_genes = load_genes(BED_PLUS, "+") + load_genes(BED_MINUS, "-")
gene_dict = {g["name"]: g for g in all_genes}
selected = selected & set(gene_dict.keys())
print(f"  After matching annotation: {len(selected)}")

# count reads for selected genes
print("  Counting reads...")
count_matrix = {}
for sample_name, bam_path in SAMPLES_ORDER:
    print(f"    {sample_name}...")
    bam = pysam.AlignmentFile(bam_path, "rb")
    counts = {}
    for gname in selected:
        g = gene_dict[gname]
        count = 0
        for read in bam.fetch(g["chrom"], g["start"], g["end"]):
            if g["strand"] == "+" and read.is_reverse:
                count += 1
            elif g["strand"] == "-" and not read.is_reverse:
                count += 1
        counts[gname] = count
    bam.close()
    count_matrix[sample_name] = counts

count_df = pd.DataFrame(count_matrix)

# normalize: log2(CPM + 1) then z-score per gene
lib_sizes = count_df.sum(axis=0)
cpm = count_df.div(lib_sizes, axis=1) * 1e6
log_cpm = np.log2(cpm + 1)
z_scores = log_cpm.subtract(log_cpm.mean(axis=1), axis=0).div(
    log_cpm.std(axis=1), axis=0
)
z_scores = z_scores.dropna()

# determine direction and cluster within groups
gene_direction = {}
for _, row in top_genes_df.iterrows():
    name = row[gene_col]
    if name in z_scores.index:
        gene_direction[name] = "up" if row["log2FoldChange"] > 0 else "down"

up_genes = [g for g in z_scores.index if gene_direction.get(g) == "up"]
down_genes = [g for g in z_scores.index if gene_direction.get(g) == "down"]


def cluster_order(gene_list, data):
    """Hierarchical clustering to order genes within a group."""
    if len(gene_list) <= 1:
        return gene_list
    sub = data.loc[gene_list]
    link = linkage(pdist(sub.values), method="ward")
    dendro = dendrogram(link, no_plot=True)
    return [gene_list[i] for i in dendro["leaves"]]


up_ordered = cluster_order(up_genes, z_scores)
down_ordered = cluster_order(down_genes, z_scores)
z_ordered = z_scores.loc[up_ordered + down_ordered]

# plot heatmap
print("  Plotting heatmap...")
n_genes = len(z_ordered)
fig_height = max(6, n_genes * 0.12 + 2)
fig, ax = plt.subplots(figsize=(6, fig_height))

cmap = sns.diverging_palette(240, 10, as_cmap=True)
im = ax.imshow(
    z_ordered.values,
    aspect="auto",
    cmap=cmap,
    vmin=-2.5,
    vmax=2.5,
    interpolation="none",
)

sample_labels = [
    "OsTir1\nCtrl 1",
    "OsTir1\nCtrl 2",
    "OsTir1\nAux 1",
    "OsTir1\nAux 2",
    "Chd1\nCtrl 1",
    "Chd1\nCtrl 2",
    "Chd1\nAux 1",
    "Chd1\nAux 2",
]
ax.set_xticks(range(8))
ax.set_xticklabels(sample_labels, fontsize=6.5, rotation=0, ha="center")
ax.xaxis.tick_top()
ax.set_yticks(range(n_genes))
ax.set_yticklabels(z_ordered.index, fontsize=4.5)
ax.yaxis.tick_right()

# color bar for conditions
cond_colors = ["#5E81AC"] * 2 + ["#81A1C1"] * 2 + ["#A3BE8C"] * 2 + ["#BF616A"] * 2
for i, color in enumerate(cond_colors):
    ax.add_patch(plt.Rectangle((i - 0.5, -1.8), 1, 1, color=color, clip_on=False))

# direction markers on left side
for i, gene in enumerate(up_ordered + down_ordered):
    d = gene_direction.get(gene, "unknown")
    color = "#BF616A" if d == "up" else "#5E81AC"
    ax.add_patch(plt.Rectangle((-0.9, i - 0.5), 0.3, 1, color=color, clip_on=False))

# separator line
if up_ordered and down_ordered:
    ax.axhline(len(up_ordered) - 0.5, color="black", linewidth=1)
    ax.text(
        -1.5,
        len(up_ordered) / 2,
        "↑",
        ha="center",
        va="center",
        fontsize=10,
        color="#BF616A",
        fontweight="bold",
        clip_on=False,
    )
    ax.text(
        -1.5,
        len(up_ordered) + len(down_ordered) / 2,
        "↓",
        ha="center",
        va="center",
        fontsize=10,
        color="#5E81AC",
        fontweight="bold",
        clip_on=False,
    )

cbar = fig.colorbar(im, ax=ax, shrink=0.4, aspect=20, pad=0.12)
cbar.set_label("Z-score", fontsize=8)
ax.set_title(f"Top {TOP_N} up + {TOP_N} down genes", pad=25)

plt.tight_layout()
fig.savefig(f"{OUT}/heatmap_top_genes.png", dpi=300, facecolor="white")
plt.close()


# --- sample correlation heatmap ---
print("  Plotting sample correlation...")
corr = log_cpm.corr()
fig, ax = plt.subplots(figsize=(5.5, 5))
im = ax.imshow(corr.values, cmap="RdYlBu_r", vmin=0.85, vmax=1.0, aspect="equal")

ax.set_xticks(range(8))
ax.set_xticklabels(sample_labels, fontsize=6.5, rotation=45, ha="right")
ax.set_yticks(range(8))
ax.set_yticklabels(sample_labels, fontsize=6.5)

for i in range(8):
    for j in range(8):
        val = corr.values[i, j]
        ax.text(
            j,
            i,
            f"{val:.3f}",
            ha="center",
            va="center",
            fontsize=5.5,
            color="white" if val < 0.93 else "black",
        )

# highlight condition blocks
for start, end, color in [
    (0, 2, "#5E81AC"),
    (2, 4, "#81A1C1"),
    (4, 6, "#A3BE8C"),
    (6, 8, "#BF616A"),
]:
    rect = plt.Rectangle(
        (start - 0.5, start - 0.5),
        end - start,
        end - start,
        fill=False,
        edgecolor=color,
        linewidth=2,
        clip_on=False,
    )
    ax.add_patch(rect)

cbar = fig.colorbar(im, ax=ax, shrink=0.7, aspect=25)
cbar.set_label("Pearson r", fontsize=8)
ax.set_title("Sample correlation", pad=10)

plt.tight_layout()
fig.savefig(f"{OUT}/sample_correlation.png", dpi=300, facecolor="white")
plt.close()

print(f"\nAll plots saved to {OUT}/")
