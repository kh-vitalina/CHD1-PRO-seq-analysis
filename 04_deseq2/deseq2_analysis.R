#!/usr/bin/env Rscript
# DESeq2 differential transcription analysis of PRO-seq data
# 4 conditions: OsTir1 ctrl/auxin, Chd1 ctrl/auxin (2 replicates each)

library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(DESeq2)
library(ggplot2)
library(pheatmap)

output_dir <- "deseq2_results"
dir.create(output_dir, showWarnings = FALSE)

# sample info
samples <- data.frame(
  sample = c("OsTir1_Ctrl_1", "OsTir1_Ctrl_2",
             "OsTir1_Auxin_1", "OsTir1_Auxin_2",
             "Chd1_Ctrl_1", "Chd1_Ctrl_2",
             "Chd1_Auxin_1", "Chd1_Auxin_2"),
  bam = c("dedup/SRR15904219_dedup.bam", "dedup/SRR15904221_dedup.bam",
           "dedup/SRR15904220_dedup.bam", "dedup/SRR15904222_dedup.bam",
           "dedup/SRR15904237_dedup.bam", "dedup/SRR15904239_dedup.bam",
           "dedup/SRR15904238_dedup.bam", "dedup/SRR15904240_dedup.bam"),
  genotype = c("OsTir1", "OsTir1", "OsTir1", "OsTir1",
               "Chd1_AID", "Chd1_AID", "Chd1_AID", "Chd1_AID"),
  treatment = c("Ctrl", "Ctrl", "Auxin", "Auxin",
                "Ctrl", "Ctrl", "Auxin", "Auxin"),
  condition = c("OsTir1_Ctrl", "OsTir1_Ctrl",
                "OsTir1_Auxin", "OsTir1_Auxin",
                "Chd1_Ctrl", "Chd1_Ctrl",
                "Chd1_Auxin", "Chd1_Auxin"),
  stringsAsFactors = FALSE
)

cat("Samples:\n")
print(samples[, c("sample", "condition")])


# load gene annotation from BED
load_bed <- function(bed_file) {
  bed <- read.table(bed_file, header = FALSE, sep = "\t",
                    quote = "", comment.char = "",
                    col.names = c("chrom", "start", "end", "name", "score", "strand"))
  bed <- bed[bed$chrom %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"), ]
  bed <- bed[(bed$end - bed$start) >= 1000, ]
  bed <- bed[bed$strand %in% c("+", "-"), ]
  return(bed)
}

genes_plus <- load_bed("annotation/genes_dm6_plus.bed")
genes_minus <- load_bed("annotation/genes_dm6_minus.bed")
all_genes <- rbind(genes_plus, genes_minus)
cat("Total genes:", nrow(all_genes), "\n")


# count reads per gene (strand-specific, PRO-seq is reverse-stranded)
count_reads <- function(bam_file, genes_bed) {
  cat("  Counting:", bam_file, "\n")
  counts <- integer(nrow(genes_bed))
  bam <- BamFile(bam_file)

  for (i in seq_len(nrow(genes_bed))) {
    chrom <- genes_bed$chrom[i]
    start <- genes_bed$start[i]
    end <- genes_bed$end[i]
    strand <- genes_bed$strand[i]

    if (strand == "+") {
      param <- ScanBamParam(
        which = GRanges(chrom, IRanges(start + 1, end)),
        flag = scanBamFlag(isMinusStrand = TRUE))
    } else {
      param <- ScanBamParam(
        which = GRanges(chrom, IRanges(start + 1, end)),
        flag = scanBamFlag(isMinusStrand = FALSE))
    }
    counts[i] <- countBam(bam, param = param)$records
  }
  return(counts)
}

# build count matrix
count_matrix <- matrix(0, nrow = nrow(all_genes), ncol = nrow(samples))
rownames(count_matrix) <- all_genes$name
colnames(count_matrix) <- samples$sample

for (j in seq_len(nrow(samples))) {
  count_matrix[, j] <- count_reads(samples$bam[j], all_genes)
  cat("    Total counts:", sum(count_matrix[, j]), "\n")
}

write.csv(count_matrix, file.path(output_dir, "raw_counts.csv"))

# filter low-expression genes
keep <- rowSums(count_matrix >= 10) >= 2
count_matrix_filtered <- count_matrix[keep, ]
cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(count_matrix_filtered), "\n")


# run DESeq2
col_data <- data.frame(
  condition = factor(samples$condition,
                     levels = c("OsTir1_Ctrl", "OsTir1_Auxin",
                                "Chd1_Ctrl", "Chd1_Auxin")),
  row.names = samples$sample
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = col_data,
  design = ~ condition
)
dds <- DESeq(dds)
cat("DESeq2 done.\n")


# PCA
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#2166AC", "#4393C3", "#1A9641", "#D6604D")) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of PRO-seq samples") +
  theme_bw(base_size = 14)
ggsave(file.path(output_dir, "PCA_plot.png"), p, width = 8, height = 6, dpi = 150)


# pairwise contrasts
contrasts <- list(
  list(name = "Auxin_effect_OsTir1",
       num = "OsTir1_Auxin", denom = "OsTir1_Ctrl"),
  list(name = "Chd1_depletion",
       num = "Chd1_Auxin", denom = "OsTir1_Auxin"),
  list(name = "Auxin_effect_Chd1",
       num = "Chd1_Auxin", denom = "Chd1_Ctrl")
)

all_results <- list()
for (comp in contrasts) {
  res <- results(dds, contrast = c("condition", comp$num, comp$denom), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$chrom <- all_genes$chrom[match(res_df$gene, all_genes$name)]
  all_results[[comp$name]] <- res_df

  sig_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
  sig_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("\n%s: up=%d, down=%d, median log2FC=%.3f\n",
              comp$name, sig_up, sig_down, median(res_df$log2FoldChange, na.rm = TRUE)))

  write.csv(res_df[order(res_df$padj), ],
            file.path(output_dir, paste0("DESeq2_", comp$name, ".csv")),
            row.names = FALSE)
}


# MA plots
png(file.path(output_dir, "MA_plots.png"), width = 1200, height = 400, res = 150)
par(mfrow = c(1, 3))
for (comp in contrasts) {
  res <- results(dds, contrast = c("condition", comp$num, comp$denom), alpha = 0.05)
  plotMA(res, main = comp$name, ylim = c(-4, 4), cex = 0.5)
}
dev.off()


# volcano plots
png(file.path(output_dir, "Volcano_plots.png"), width = 1200, height = 400, res = 150)
par(mfrow = c(1, 3))
for (comp in contrasts) {
  res_df <- all_results[[comp$name]]
  res_df <- res_df[!is.na(res_df$padj), ]
  cols <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "red",
          ifelse(res_df$padj < 0.05, "orange", "grey"))
  plot(res_df$log2FoldChange, -log10(res_df$padj),
       pch = 16, cex = 0.4, col = cols,
       xlab = "log2FC", ylab = "-log10(padj)",
       main = comp$name, xlim = c(-5, 5))
  abline(h = -log10(0.05), col = "blue", lty = 2)
  abline(v = c(-1, 1), col = "blue", lty = 2)
}
dev.off()


# sample correlation heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
annotation_col <- data.frame(Condition = col_data$condition, row.names = rownames(col_data))
ann_colors <- list(
  Condition = c(OsTir1_Ctrl = "#2166AC", OsTir1_Auxin = "#4393C3",
                Chd1_Ctrl = "#1A9641", Chd1_Auxin = "#D6604D")
)

png(file.path(output_dir, "sample_correlation_heatmap.png"), width = 800, height = 700, res = 150)
pheatmap(sample_dist_matrix, annotation_col = annotation_col,
         annotation_colors = ann_colors, main = "Sample distance matrix")
dev.off()


# X vs autosomes for main comparison
res_main <- all_results[["Chd1_depletion"]]
res_main <- res_main[!is.na(res_main$padj), ]
chrX <- res_main[res_main$chrom == "chrX", ]
auto <- res_main[res_main$chrom != "chrX", ]

cat("\nChd1 depletion — X vs autosomes:\n")
cat(sprintf("  X (%d genes): median log2FC = %.3f, sig up=%d, sig down=%d\n",
    nrow(chrX), median(chrX$log2FoldChange, na.rm = TRUE),
    sum(chrX$padj < 0.05 & chrX$log2FoldChange > 0, na.rm = TRUE),
    sum(chrX$padj < 0.05 & chrX$log2FoldChange < 0, na.rm = TRUE)))
cat(sprintf("  Autosomes (%d genes): median log2FC = %.3f, sig up=%d, sig down=%d\n",
    nrow(auto), median(auto$log2FoldChange, na.rm = TRUE),
    sum(auto$padj < 0.05 & auto$log2FoldChange > 0, na.rm = TRUE),
    sum(auto$padj < 0.05 & auto$log2FoldChange < 0, na.rm = TRUE)))


# summary
cat("\nSummary:\n")
for (comp in contrasts) {
  res_df <- all_results[[comp$name]]
  sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  cat(sprintf("  %s: %d significant genes\n", comp$name, sig))
}
cat("\nResults saved to", output_dir, "\n")
