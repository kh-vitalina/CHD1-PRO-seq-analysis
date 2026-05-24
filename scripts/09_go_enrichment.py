#!/usr/bin/env python3
"""
GO enrichment analysis for DESeq2 results
Sends up/downregulated gene lists to g:Profiler API,
plots dot plots and combined barplot
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import requests
import os

OUT = 'deseq2_results'
os.makedirs(OUT, exist_ok=True)


def query_gprofiler(gene_list, label='query'):
    """Send gene list to g:Profiler REST API, return enriched terms."""
    print(f"  Querying g:Profiler for {label} ({len(gene_list)} genes)...")
    url = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile/'
    payload = {
        'organism': 'dmelanogaster',
        'query': list(gene_list),
        'sources': ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG'],
        'user_threshold': 0.05,
        'significance_threshold_method': 'fdr',
        'no_evidences': True,
    }
    try:
        resp = requests.post(url, json=payload, timeout=60)
        resp.raise_for_status()
        results = resp.json().get('result', [])
    except Exception as e:
        print(f"  Request failed: {e}")
        return pd.DataFrame()

    if not results:
        print(f"  No significant terms for {label}")
        return pd.DataFrame()

    rows = [{
        'source': r['source'],
        'term': r['name'],
        'term_id': r['native'],
        'padj': r['p_value'],
        'term_size': r['term_size'],
        'intersection_size': r['intersection_size'],
        'query_size': r['query_size'],
        'ratio': r['intersection_size'] / r['query_size'],
    } for r in results]

    df = pd.DataFrame(rows)
    df['-log10(padj)'] = -np.log10(df['padj'])
    print(f"  Found {len(df)} significant terms")
    return df


def plot_dotplot(go_df, title, filename, top_n=15):
    """Dot plot: y = GO terms, x = gene ratio, size = gene count, color = source."""
    if len(go_df) == 0:
        print(f"  Skipping {filename}, no data")
        return

    source_colors = {
        'GO:BP': '#E74C3C', 'GO:MF': '#3498DB',
        'GO:CC': '#2ECC71', 'KEGG': '#9B59B6',
    }

    top = go_df.sort_values('padj').drop_duplicates('term').head(top_n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6, max(3, len(top) * 0.35 + 1)))

    colors = [source_colors.get(s, '#888') for s in top['source']]
    ax.scatter(top['ratio'], range(len(top)), s=top['intersection_size'] * 12,
               c=colors, alpha=0.75, edgecolors='white', linewidth=0.8)

    labels = [t[:55] + '...' if len(t) > 55 else t for t in top['term']]
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('Gene ratio (intersection / query)')
    ax.set_title(title)

    # source legend
    handles = [mpatches.Patch(color=c, label=s) for s, c in source_colors.items()
               if s in top['source'].values]
    if handles:
        ax.legend(handles=handles, loc='lower right', fontsize=7)

    plt.tight_layout()
    fig.savefig(f'{OUT}/{filename}.png', dpi=300, facecolor='white')
    plt.close()
    print(f"  Saved: {filename}.png")


# load DESeq2 results
deseq2_file = 'deseq2_results/DESeq2_Chd1_depletion.csv'
print("Loading DESeq2 results...")
if not os.path.exists(deseq2_file):
    print(f"ERROR: {deseq2_file} not found, run deseq2_analysis.R first")
    exit(1)

df = pd.read_csv(deseq2_file).dropna(subset=['log2FoldChange', 'padj'])

# find gene column
gene_col = 'gene'
if gene_col not in df.columns:
    for col in ['gene_name', 'name', 'Gene']:
        if col in df.columns:
            gene_col = col
            break
    else:
        df['gene'] = df.index
        gene_col = 'gene'

up_genes = df[(df['padj'] < 0.05) & (df['log2FoldChange'] >= 1)]
down_genes = df[(df['padj'] < 0.05) & (df['log2FoldChange'] <= -1)]
print(f"  Up: {len(up_genes)}, Down: {len(down_genes)}")

# run enrichment
go_up = query_gprofiler(up_genes[gene_col].tolist(), label='upregulated')
go_down = query_gprofiler(down_genes[gene_col].tolist(), label='downregulated')

if len(go_up) > 0:
    go_up.to_csv(f'{OUT}/GO_upregulated.csv', index=False)
if len(go_down) > 0:
    go_down.to_csv(f'{OUT}/GO_downregulated.csv', index=False)

# dot plots
print("\nPlotting...")
plot_dotplot(go_up, 'GO — Upregulated upon Chd1 depletion', 'GO_upregulated')
plot_dotplot(go_down, 'GO — Downregulated upon Chd1 depletion', 'GO_downregulated')

# combined barplot (top 10 each)
if len(go_up) > 0 and len(go_down) > 0:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for ax, go_df, title, color in [
        (ax1, go_up, 'Upregulated', '#BF616A'),
        (ax2, go_down, 'Downregulated', '#5E81AC'),
    ]:
        top10 = go_df.sort_values('padj').drop_duplicates('term').head(10).iloc[::-1]
        ax.barh(range(len(top10)), top10['-log10(padj)'], color=color,
                alpha=0.75, edgecolor='white', height=0.65)
        for i, (_, row) in enumerate(top10.iterrows()):
            ax.text(row['-log10(padj)'] + 0.1, i,
                    f"{row['intersection_size']} genes", va='center', fontsize=6.5, color='#555')
        labels = [t[:45] + '...' if len(t) > 45 else t for t in top10['term']]
        ax.set_yticks(range(len(top10)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel(r'$-\log_{10}$(adjusted p-value)')
        ax.set_title(title, color=color)

    plt.tight_layout()
    fig.savefig(f'{OUT}/GO_combined.png', dpi=300, facecolor='white')
    plt.close()
    print(f"  Saved: GO_combined.png")
    