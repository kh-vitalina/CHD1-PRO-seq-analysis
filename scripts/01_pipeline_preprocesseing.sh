#!/bin/bash
# PRO-seq preprocessing pipeline: FASTQ -> BAM -> BigWig -> metagene
# Drosophila S2 cells, dm6 genome
# Data: Hendy et al. 2022 (GEO: GSE184187)

set -euo pipefail

GENOME_DIR="genome"
GENOME_FA="$GENOME_DIR/dm6.fa"
INDEX="$GENOME_DIR/dm6_index"
THREADS=4

SAMPLES=(
    SRR15904219   # OsTir1_Ctrl_1
    SRR15904221   # OsTir1_Ctrl_2
    SRR15904220   # OsTir1_Auxin_1
    SRR15904222   # OsTir1_Auxin_2
    SRR15904237   # Chd1_Ctrl_1
    SRR15904239   # Chd1_Ctrl_2
    SRR15904238   # Chd1_Auxin_1
    SRR15904240   # Chd1_Auxin_2
)

mkdir -p fastq trimmed_adapters aligned dedup bigwig annotation qc/raw qc/trimmed


# 1. Download from SRA
echo "--- Downloading FASTQ ---"
for s in "${SAMPLES[@]}"; do
    if [ ! -f "fastq/${s}.fastq.gz" ]; then
        fasterq-dump --split-files -O fastq "$s"
        gzip "fastq/${s}.fastq"
    fi
done


# 2. FastQC on raw reads
echo "--- FastQC (raw) ---"
fastqc -o qc/raw -t $THREADS fastq/*.fastq.gz


# 3. UMI extraction (8 nt barcode at 5' end)
echo "--- UMI extraction ---"
for s in "${SAMPLES[@]}"; do
    if [ ! -f "fastq/${s}_umi.fastq.gz" ]; then
        umi_tools extract \
            --bc-pattern=NNNNNNNN \
            -I "fastq/${s}.fastq.gz" \
            -S "fastq/${s}_umi.fastq.gz"
    fi
done


# 4. Adapter trimming
echo "--- Trimming adapters ---"
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
for s in "${SAMPLES[@]}"; do
    if [ ! -f "trimmed_adapters/${s}_trimmed.fastq.gz" ]; then
        cutadapt \
            -a "$ADAPTER" \
            -m 15 -O 10 \
            -o "trimmed_adapters/${s}_trimmed.fastq.gz" \
            "fastq/${s}_umi.fastq.gz" \
            > "trimmed_adapters/${s}_cutadapt.log"
    fi
done


# 5. FastQC on trimmed reads
echo "--- FastQC (trimmed) ---"
fastqc -o qc/trimmed -t $THREADS trimmed_adapters/*_trimmed.fastq.gz


# 6. Build bowtie2 index (if needed)
if [ ! -f "${INDEX}.1.bt2" ]; then
    echo "--- Building bowtie2 index ---"
    bowtie2-build "$GENOME_FA" "$INDEX"
fi


# 7. Alignment with bowtie2
echo "--- Aligning with bowtie2 ---"
for s in "${SAMPLES[@]}"; do
    if [ ! -f "aligned/${s}.bam" ]; then
        bowtie2 \
            -x "$INDEX" \
            -U "trimmed_adapters/${s}_trimmed.fastq.gz" \
            -p $THREADS \
            --very-sensitive \
            2> "aligned/${s}_bowtie2.log" \
        | samtools sort -@ $THREADS -o "aligned/${s}.bam"
        samtools index "aligned/${s}.bam"
    fi
done


# 8. UMI deduplication
echo "--- UMI deduplication ---"
for s in "${SAMPLES[@]}"; do
    if [ ! -f "dedup/${s}_dedup.bam" ]; then
        umi_tools dedup \
            -I "aligned/${s}.bam" \
            -S "dedup/${s}_dedup.bam" \
            --output-stats="dedup/${s}_stats" \
            > "dedup/${s}_dedup.log" 2>&1
        samtools index "dedup/${s}_dedup.bam"
    fi
done


# 9. Strand-specific BigWig generation (bedtools method)
# PRO-seq: reverse reads = plus strand signal, forward reads = minus strand
echo "--- Generating BigWig files ---"
CHROM_SIZES="$GENOME_DIR/dm6.chrom.sizes"
if [ ! -f "$CHROM_SIZES" ]; then
    samtools faidx "$GENOME_FA"
    cut -f1,2 "${GENOME_FA}.fai" > "$CHROM_SIZES"
fi

mkdir -p bigwig
for s in "${SAMPLES[@]}"; do
    BAM="dedup/${s}_dedup.bam"

    # plus strand signal (from reverse reads, shift +1)
    if [ ! -f "bigwig/${s}_ps.bw" ]; then
        samtools view -b -f 16 "$BAM" \
        | bedtools genomecov -bg -ibam stdin \
        | awk 'BEGIN{OFS="\t"}{print $1, $2+1, $3+1, $4}' \
        | sort -k1,1 -k2,2n \
        | bedGraphToBigWig stdin "$CHROM_SIZES" "bigwig/${s}_ps.bw"
    fi

    # minus strand signal (from forward reads, shift -1)
    if [ ! -f "bigwig/${s}_ns.bw" ]; then
        samtools view -b -F 16 "$BAM" \
        | bedtools genomecov -bg -ibam stdin \
        | awk 'BEGIN{OFS="\t"}{if($2>0) print $1, $2-1, $3-1, $4}' \
        | sort -k1,1 -k2,2n \
        | bedGraphToBigWig stdin "$CHROM_SIZES" "bigwig/${s}_ns.bw"
    fi
done


# 10. Average replicates (BigWig)
echo "--- Averaging replicates ---"

avg_bigwig() {
    local out=$1; shift
    local bw1=$1; shift
    local bw2=$1

    if [ ! -f "$out" ]; then
        bigwigMerge "$bw1" "$bw2" /dev/stdout \
        | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4/2}' \
        | sort -k1,1 -k2,2n \
        | bedGraphToBigWig stdin "$CHROM_SIZES" "$out"
    fi
}

avg_bigwig bigwig/OsTir1_Ctrl_mean_ps.bw   bigwig/SRR15904219_ps.bw bigwig/SRR15904221_ps.bw
avg_bigwig bigwig/OsTir1_Auxin_mean_ps.bw  bigwig/SRR15904220_ps.bw bigwig/SRR15904222_ps.bw
avg_bigwig bigwig/Chd1_Ctrl_mean_ps.bw     bigwig/SRR15904237_ps.bw bigwig/SRR15904239_ps.bw
avg_bigwig bigwig/Chd1_Auxin_mean_ps.bw    bigwig/SRR15904238_ps.bw bigwig/SRR15904240_ps.bw

avg_bigwig bigwig/OsTir1_Ctrl_mean_ns.bw   bigwig/SRR15904219_ns.bw bigwig/SRR15904221_ns.bw
avg_bigwig bigwig/OsTir1_Auxin_mean_ns.bw  bigwig/SRR15904220_ns.bw bigwig/SRR15904222_ns.bw
avg_bigwig bigwig/Chd1_Ctrl_mean_ns.bw     bigwig/SRR15904237_ns.bw bigwig/SRR15904239_ns.bw
avg_bigwig bigwig/Chd1_Auxin_mean_ns.bw    bigwig/SRR15904238_ns.bw bigwig/SRR15904240_ns.bw


# 11. Metagene matrix (deepTools)
echo "--- Computing metagene matrix ---"

# plus strand genes, plus strand signal
computeMatrix scale-regions \
    -S bigwig/OsTir1_Ctrl_mean_ps.bw \
       bigwig/OsTir1_Auxin_mean_ps.bw \
       bigwig/Chd1_Ctrl_mean_ps.bw \
       bigwig/Chd1_Auxin_mean_ps.bw \
    -R annotation/genes_dm6_plus.bed \
    -b 500 -a 500 --regionBodyLength 2000 \
    --binSize 10 --missingDataAsZero \
    -o matrix_plus.gz -p $THREADS

# minus strand genes, minus strand signal
computeMatrix scale-regions \
    -S bigwig/OsTir1_Ctrl_mean_ns.bw \
       bigwig/OsTir1_Auxin_mean_ns.bw \
       bigwig/Chd1_Ctrl_mean_ns.bw \
       bigwig/Chd1_Auxin_mean_ns.bw \
    -R annotation/genes_dm6_minus.bed \
    -b 500 -a 500 --regionBodyLength 2000 \
    --binSize 10 --missingDataAsZero \
    -o matrix_minus.gz -p $THREADS

echo "--- Done ---"
echo "Next steps:"
echo "  1. Run scripts/02_prepare_annotation.sh"
echo "  2. Run scripts/03_plot_metagene.py"
echo "  3. Run analysis scripts (04-10)"
