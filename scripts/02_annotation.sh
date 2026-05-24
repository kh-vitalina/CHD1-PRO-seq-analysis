#!/bin/bash
# Prepare gene annotation from UCSC refGene for PRO-seq analysis
# Downloads refGene table, converts to BED, splits by strand

set -euo pipefail

mkdir -p annotation

REFGENE="annotation/refGene_dm6.txt"
OUTBED="annotation/genes_dm6.bed"

# download refGene if needed
if [ ! -f "$REFGENE" ]; then
    echo "Downloading refGene from UCSC..."
    curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/refGene.txt.gz" \
    | gunzip > "$REFGENE"
fi

# convert to BED6 (chrom, start, end, name, score, strand)
# refGene columns: bin, name, chrom, strand, txStart, txEnd, ...
# take longest transcript per gene
echo "Converting to BED..."
awk -F'\t' 'BEGIN{OFS="\t"} {
    gene = $13
    chrom = $3
    strand = $4
    start = $5
    end = $6
    len = end - start
    key = gene
    if (!(key in best) || len > best_len[key]) {
        best[key] = chrom"\t"start"\t"end"\t"gene"\t0\t"strand
        best_len[key] = len
    }
}
END { for (k in best) print best[k] }' "$REFGENE" \
| sort -k1,1 -k2,2n > "$OUTBED"

# split by strand
awk '$6 == "+"' "$OUTBED" > annotation/genes_dm6_plus.bed
awk '$6 == "-"' "$OUTBED" > annotation/genes_dm6_minus.bed

# unique genes only (no duplicates by name)
sort -k4,4 -u "$OUTBED" > annotation/genes_dm6_unique.bed

echo "Done:"
echo "  $(wc -l < annotation/genes_dm6_plus.bed) plus-strand genes"
echo "  $(wc -l < annotation/genes_dm6_minus.bed) minus-strand genes"
echo "  $(wc -l < annotation/genes_dm6_unique.bed) unique genes"
