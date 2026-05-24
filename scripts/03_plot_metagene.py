#!/usr/bin/env python3
"""
Combined metagene plot from deepTools matrices
Averages plus and minus strand profiles weighted by gene count
"""

import numpy as np
import matplotlib.pyplot as plt
import gzip
import json


def read_deeptools_matrix(matrix_file):
    """Parse deepTools computeMatrix output, return mean profile per sample"""
    with gzip.open(matrix_file, 'rt') as f:
        header = json.loads(f.readline().strip().replace('@', ''))
        sample_labels = header['sample_labels']
        boundaries = header['sample_boundaries']
        data = []
        for line in f:
            vals = line.strip().split('\t')
            data.append([float(x) if x != 'nan' else np.nan for x in vals[6:]])

    data = np.array(data)
    profiles = {}
    for i, label in enumerate(sample_labels):
        profiles[label] = np.nanmean(data[:, boundaries[i]:boundaries[i+1]], axis=0)
    return profiles, len(data)


# read both strand matrices
plus_profiles, n_plus = read_deeptools_matrix('matrix_bowtie1_plus_all.gz')
minus_profiles, n_minus = read_deeptools_matrix('matrix_bowtie1_minus_pos.gz')
print(f"Plus strand genes: {n_plus}, Minus strand genes: {n_minus}")

# sample pairs: (plus key, minus key, display name, color)
sample_pairs = [
    ('OsTir1_Ctrl_mean_ps',  'OsTir1_Ctrl_mean_ns_pos',  'OsTir1_Ctrl',  '#2166AC'),
    ('OsTir1_Auxin_mean_ps', 'OsTir1_Auxin_mean_ns_pos', 'OsTir1_Auxin', '#D6604D'),
    ('Chd1_Ctrl_mean_ps',    'Chd1_Ctrl_mean_ns_pos',    'Chd1_Ctrl',    '#1A9641'),
    ('Chd1_Auxin_mean_ps',   'Chd1_Auxin_mean_ns_pos',   'Chd1_Auxin',   '#7B2D8B'),
]

fig, ax = plt.subplots(figsize=(8, 5))

for ps_key, ns_key, label, color in sample_pairs:
    p = plus_profiles[ps_key]
    m = minus_profiles[ns_key]
    n = min(len(p), len(m))
    # weighted average by number of genes on each strand
    combined = (p[:n] * n_plus + m[:n] * n_minus) / (n_plus + n_minus)
    x = np.linspace(-0.5, 2.5, len(combined))
    ax.plot(x, combined, color=color, label=label, linewidth=1.5)

ax.axvline(0, color='black', ls='--', alpha=0.3)
ax.axvline(2.0, color='black', ls='--', alpha=0.3)
ax.text(0, ax.get_ylim()[1] * 0.95, 'TSS', ha='center', fontsize=10)
ax.text(2.0, ax.get_ylim()[1] * 0.95, 'TES', ha='center', fontsize=10)
ax.set_xlabel('Distance (kb)')
ax.set_ylabel('PRO-seq signal (TPM)')
ax.set_title('PRO-seq metagene: all genes (bowtie, dm6)')
ax.legend()
plt.tight_layout()
plt.savefig('metagene_combined.png', dpi=150)
plt.close()
print("Saved: metagene_combined.png")
