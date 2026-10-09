#!/usr/bin/env Rscript

# STAR + featureCounts differential-expression analysis for A. nidulans.
# Usage:
# Rscript 06_deseq2_analysis.R metadata/sample_key.csv results/counts/gene_counts.txt

required_packages <- c(
  "DESeq2", "apeglm", "dplyr", "tidyr", "tibble", "ggplot2",
  "purrr", "readr", "pheatmap", "vsn", "RColorBrewer", "ggrepel"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(pheatmap)
  library(vsn)
  library(RColorBrewer)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "Usage: Rscript 06_deseq2_analysis.R metadata/sample_key.csv results/counts/gene_counts.txt",
    call. = FALSE
  )
}
metadata_file <- args[1]
count_matrix_file <- args[2]

alpha <- 0.05
log2fc_threshold <- 1
min_count <- 10
min_samples <- 3
wt_reference <- "WT"
output_dir <- "results/deseq2"

custom_comparisons <- tribble(
  ~numerator, ~denominator,
  "double_del", "leuB_del",
  "double_del", "leuR_del"
)

for (input_file in c(metadata_file, count_matrix_file)) {
  if (!file.exists(input_file)) stop("Input not found: ", input_file, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
qc_dir <- file.path(output_dir, "qc")
results_dir <- file.path(output_dir, "results")
significant_dir <- file.path(output_dir, "significant")
plots_dir <- file.path(output_dir, "plots")
ma_dir <- file.path(plots_dir, "MA")
volcano_dir <- file.path(plots_dir, "volcano")
for (path in c(qc_dir, results_dir, significant_dir, ma_dir, volcano_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

message("Inputs validated successfully.")
message("Loading count matrix...")
counts_raw <- read.delim(
  count_matrix_file,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_fc_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
if (!all(required_fc_columns %in% names(counts_raw))) {
  stop("Input is not a valid featureCounts matrix.", call. = FALSE)
}
if (anyDuplicated(counts_raw$Geneid)) stop("Duplicate Geneid values in count matrix.", call. = FALSE)

counts <- counts_raw[, setdiff(names(counts_raw), required_fc_columns[-1]), drop = FALSE]
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL
counts <- as.matrix(counts)
if (anyNA(counts)) stop("Count matrix contains missing values.", call. = FALSE)
if (any(counts < 0)) stop("Count matrix contains negative values.", call. = FALSE)
storage.mode(counts) <- "integer"

clean_sample_names <- function(x) {
  x <- basename(x)
  sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", x)
}
colnames(counts) <- clean_sample_names(colnames(counts))
if (anyDuplicated(colnames(counts))) stop("Duplicate cleaned sample names in count matrix.", call. = FALSE)
message("Loaded ", nrow(counts), " genes across ", ncol(counts), " samples.")

message("Loading metadata...")
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
required_metadata_columns <- c("sample", "genotype", "nitrogen", "replicate", "batch")
missing_columns <- setdiff(required_metadata_columns, names(metadata))
if (length(missing_columns)) {
  stop("Metadata missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
metadata <- metadata[, required_metadata_columns, drop = FALSE]
if (anyDuplicated(metadata$sample)) stop("Duplicate sample names in metadata.", call. = FALSE)

count_samples <- colnames(counts)
metadata_samples <- metadata$sample
if (length(setdiff(count_samples, metadata_samples))) {
  stop("Count samples absent from metadata: ", paste(setdiff(count_samples, metadata_samples), collapse = ", "), call. = FALSE)
}
if (length(setdiff(metadata_samples, count_samples))) {
  stop("Metadata samples absent from counts: ", paste(setdiff(metadata_samples, count_samples), collapse = ", "), call. = FALSE)
}
message("Metadata and count matrix samples match.")

metadata <- metadata[match(count_samples, metadata$sample), , drop = FALSE]
rownames(metadata) <- metadata$sample
stopifnot(identical(rownames(metadata), colnames(counts)))
metadata$genotype <- factor(metadata$genotype)
metadata$nitrogen <- factor(metadata$nitrogen)
metadata$batch <- factor(metadata$batch)
if (!wt_reference %in% levels(metadata$genotype)) stop("WT reference not found.", call. = FALSE)
metadata$genotype <- relevel(metadata$genotype, ref = wt_reference)

design_terms <- character(0)
if (nlevels(droplevels(metadata$batch)) > 1) design_terms <- c(design_terms, "batch")
if (nlevels(droplevels(metadata$nitrogen)) > 1) design_terms <- c(design_terms, "nitrogen")
design_terms <- c(design_terms, "genotype")
design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
message("Design formula: ", deparse(design_formula))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = design_formula)
genes_before_filter <- nrow(dds)
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]
if (!nrow(dds)) stop("No genes passed the count filter.", call. = FALSE)
message("Retained ", nrow(dds), " of ", genes_before_filter, " genes after filtering.")
write.csv(
  data.frame(metric = c("genes_before_filter", "genes_after_filter"), value = c(genes_before_filter, nrow(dds))),
  file.path(output_dir, "filter_summary.csv"), row.names = FALSE
)
message("DESeq2 dataset constructed successfully.")

message("Running DESeq2 model fitting...")
dds <- DESeq(dds, parallel = FALSE)
message("DESeq2 model completed.")
write.csv(data.frame(coefficient = resultsNames(dds)), file.path(output_dir, "results_names.csv"), row.names = FALSE)

export_matrix <- function(matrix_object, filename) {
  as.data.frame(matrix_object) |>
    rownames_to_column("gene_id") |>
    write.csv(file.path(output_dir, filename), row.names = FALSE)
}
export_matrix(counts(dds, normalized = FALSE), "raw_counts.csv")
message("Raw counts exported.")
export_matrix(counts(dds, normalized = TRUE), "normalized_counts.csv")
message("Normalized counts exported.")

message("Generating VST matrix...")
vsd <- vst(dds, blind = TRUE)
export_matrix(assay(vsd), "vst_counts.csv")
message("VST counts exported.")

write.csv(
  data.frame(sample = colnames(dds), size_factor = sizeFactors(dds)),
  file.path(output_dir, "size_factors.csv"), row.names = FALSE
)
write.csv(
  data.frame(gene_id = rownames(dds), dispersion = dispersions(dds)),
  file.path(output_dir, "dispersion_estimates.csv"), row.names = FALSE
)
write.csv(
  data.frame(
    metric = c("samples", "genes_after_filter", "genotypes", "nitrogen_conditions", "batches"),
    value = c(ncol(dds), nrow(dds), nlevels(metadata$genotype), nlevels(metadata$nitrogen), nlevels(metadata$batch))
  ),
  file.path(output_dir, "analysis_summary.csv"), row.names = FALSE
)
pdf(file.path(output_dir, "dispersion_plot.pdf"), width = 7, height = 6)
plotDispEsts(dds)
dev.off()
saveRDS(dds, file.path(output_dir, "dds.rds"))
saveRDS(vsd, file.path(output_dir, "vsd.rds"))
message("DESeq2 objects saved.")

message("Generating PCA plot...")
pca_groups <- "genotype"
if (nlevels(metadata$nitrogen) > 1) pca_groups <- c(pca_groups, "nitrogen")
pca_data <- plotPCA(vsd, intgroup = pca_groups, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
if (nlevels(metadata$nitrogen) > 1) {
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = genotype, shape = nitrogen))
} else {
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = genotype))
}
pca_plot <- pca_plot +
  geom_point(size = 4) +
  labs(
    title = "PCA of RNA-seq Samples",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(qc_dir, "qc_pca.pdf"), pca_plot, width = 8, height = 6)
write.csv(pca_data, file.path(qc_dir, "pca_coordinates.csv"), row.names = FALSE)
write.csv(
  data.frame(principal_component = c("PC1", "PC2"), variance_percent = percent_var),
  file.path(qc_dir, "pca_variance_explained.csv"), row.names = FALSE
)

annotation_df <- metadata[, c("genotype", "nitrogen", "batch"), drop = FALSE]
sample_distances <- dist(t(assay(vsd)))
distance_matrix <- as.matrix(sample_distances)
write.csv(distance_matrix, file.path(qc_dir, "sample_distance_matrix.csv"))
pdf(file.path(qc_dir, "qc_sample_distance_heatmap.pdf"), width = 8, height = 7)
pheatmap(
  distance_matrix,
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  clustering_distance_rows = sample_distances,
  clustering_distance_cols = sample_distances,
  main = "Sample Distance Heatmap"
)
dev.off()

correlation_matrix <- cor(assay(vsd), method = "pearson")
write.csv(correlation_matrix, file.path(qc_dir, "sample_correlation_matrix.csv"))
pdf(file.path(qc_dir, "qc_sample_correlation_heatmap.pdf"), width = 8, height = 7)
pheatmap(correlation_matrix, annotation_col = annotation_df, annotation_row = annotation_df, main = "Sample Correlation Heatmap")
dev.off()

pdf(file.path(qc_dir, "qc_sample_clustering.pdf"), width = 9, height = 7)
plot(hclust(sample_distances), main = "Hierarchical Clustering of Samples", xlab = "", sub = "")
dev.off()

sample_qc_summary <- metadata |>
  mutate(size_factor = sizeFactors(dds))
write.csv(sample_qc_summary, file.path(qc_dir, "sample_qc_summary.csv"), row.names = FALSE)
mean_correlations <- colMeans(correlation_matrix)
write.csv(
  data.frame(sample = names(mean_correlations), mean_correlation = mean_correlations) |>
    arrange(mean_correlation),
  file.path(qc_dir, "sample_mean_correlations.csv"), row.names = FALSE
)
message("QC figure generation complete.")

comparison_summary <- list()
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

classify_gene <- function(
  shrunken_log2fc,
  adjusted_pvalue,
  lfc_cutoff = log2fc_threshold,
  fdr_cutoff = alpha
) {
  dplyr::case_when(
    !is.na(adjusted_pvalue) & adjusted_pvalue < fdr_cutoff & shrunken_log2fc >= lfc_cutoff ~ "UP",
    !is.na(adjusted_pvalue) & adjusted_pvalue < fdr_cutoff & shrunken_log2fc <= -lfc_cutoff ~ "DOWN",
    TRUE ~ "NOT_SIGNIFICANT"
  )
}

create_result_table <- function(raw_result, shrunken_result, shrinkage_method) {
  raw_df <- as.data.frame(raw_result)
  shrink_df <- as.data.frame(shrunken_result)
  result_table <- tibble(
    gene_id = rownames(raw_df),
    baseMean = raw_df$baseMean,
    raw_log2FC = raw_df$log2FoldChange,
    shrunken_log2FC = shrink_df$log2FoldChange,
    lfc_difference = shrink_df$log2FoldChange - raw_df$log2FoldChange,
    shrinkage_method = shrinkage_method,
    lfcSE = shrink_df$lfcSE,
    stat = raw_df$stat,
    pvalue = raw_df$pvalue,
    padj = raw_df$padj,
    neg_log10_fdr = -log10(pmax(raw_df$padj, .Machine$double.xmin))
  )
  result_table$classification <- classify_gene(result_table$shrunken_log2FC, result_table$padj)
  arrange(result_table, padj)
}

make_ma_plot <- function(shrunken_result, comparison_label, output_file) {
  pdf(output_file, width = 7, height = 6)
  plotMA(shrunken_result, ylim = c(-5, 5), main = comparison_label)
  dev.off()
}

make_volcano_plot <- function(result_table, comparison_label, output_file) {
  label_genes <- bind_rows(
    result_table |> filter(classification == "UP") |> arrange(padj) |> slice_head(n = 15),
    result_table |> filter(classification == "DOWN") |> arrange(padj) |> slice_head(n = 15)
  )
  volcano <- ggplot(result_table, aes(shrunken_log2FC, neg_log10_fdr, color = classification)) +
    geom_point(alpha = 0.6, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey50") +
    geom_text_repel(
      data = label_genes,
      aes(label = gene_id),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      seed = 42,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(DOWN = "#2878B5", NOT_SIGNIFICANT = "grey70", UP = "#C82423")) +
    labs(title = comparison_label, x = "Shrunken log2 Fold Change", y = "-log10(FDR)", color = NULL) +
    theme_classic()
  ggsave(output_file, volcano, width = 7, height = 6)
}

write_comparison <- function(raw_result, shrunken_result, comparison_name, shrinkage_method) {
  safe_comparison <- safe_name(comparison_name)
  results_table <- create_result_table(raw_result, shrunken_result, shrinkage_method)
  significant_table <- results_table |> filter(classification != "NOT_SIGNIFICANT")
  top25_up <- results_table |>
    filter(classification == "UP") |>
    arrange(padj, desc(shrunken_log2FC)) |>
    slice_head(n = 25) |>
    mutate(rank = row_number()) |>
    relocate(rank)
  top25_down <- results_table |>
    filter(classification == "DOWN") |>
    arrange(padj, shrunken_log2FC) |>
    slice_head(n = 25) |>
    mutate(rank = row_number()) |>
    relocate(rank)
  write.csv(results_table, file.path(results_dir, paste0(safe_comparison, ".csv")), row.names = FALSE)
  write.csv(significant_table, file.path(significant_dir, paste0("significant_", safe_comparison, ".csv")), row.names = FALSE)
  write.csv(top25_up, file.path(significant_dir, paste0("top25_up_", safe_comparison, ".csv")), row.names = FALSE)
  write.csv(top25_down, file.path(significant_dir, paste0("top25_down_", safe_comparison, ".csv")), row.names = FALSE)
  comparison_summary[[comparison_name]] <<- tibble(
    comparison = comparison_name,
    total_significant = nrow(significant_table),
    upregulated = sum(significant_table$classification == "UP"),
    downregulated = sum(significant_table$classification == "DOWN")
  )
  make_ma_plot(shrunken_result, comparison_name, file.path(ma_dir, paste0("MA_", safe_comparison, ".pdf")))
  make_volcano_plot(results_table, comparison_name, file.path(volcano_dir, paste0("volcano_", safe_comparison, ".pdf")))
  message("Completed comparison: ", comparison_name)
}

message("WT COMPARISONS")
mutant_genotypes <- setdiff(levels(metadata$genotype), wt_reference)
if (!length(mutant_genotypes)) stop("No mutant genotypes found.", call. = FALSE)
for (mutant in mutant_genotypes) {
  coefficient_name <- paste0("genotype_", make.names(mutant), "_vs_", make.names(wt_reference))
  if (!coefficient_name %in% resultsNames(dds)) {
    warning("Skipping coefficient not found: ", coefficient_name)
    next
  }
  raw_result <- results(dds, name = coefficient_name, alpha = alpha)
  shrunken_result <- lfcShrink(dds, coef = coefficient_name, type = "apeglm")
  write_comparison(raw_result, shrunken_result, paste(mutant, "vs", wt_reference), "apeglm")
}

message("CUSTOM COMPARISONS")
for (i in seq_len(nrow(custom_comparisons))) {
  numerator <- custom_comparisons$numerator[i]
  denominator <- custom_comparisons$denominator[i]
  if (!all(c(numerator, denominator) %in% levels(metadata$genotype))) {
    warning("Skipping absent comparison: ", numerator, " vs ", denominator)
    next
  }
  raw_result <- results(dds, contrast = c("genotype", numerator, denominator), alpha = alpha)
  shrunken_result <- lfcShrink(
    dds,
    contrast = c("genotype", numerator, denominator),
    type = "normal"
  )
  write_comparison(raw_result, shrunken_result, paste(numerator, "vs", denominator), "normal")
}

comparison_summary_df <- bind_rows(comparison_summary)
if (!nrow(comparison_summary_df)) stop("No differential-expression comparisons completed.", call. = FALSE)
comparison_summary_df <- arrange(comparison_summary_df, desc(total_significant))
write.csv(comparison_summary_df, file.path(output_dir, "comparison_summary.csv"), row.names = FALSE)
write.csv(
  tibble(metric = c("comparisons", "samples", "genes_tested"), value = c(nrow(comparison_summary_df), ncol(dds), nrow(dds))),
  file.path(output_dir, "overall_summary.csv"), row.names = FALSE
)

comparison_summary_long <- comparison_summary_df |>
  select(comparison, upregulated, downregulated) |>
  pivot_longer(c(upregulated, downregulated), names_to = "direction", values_to = "gene_count")
deg_barplot <- ggplot(comparison_summary_long, aes(comparison, gene_count, fill = direction)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(upregulated = "#C82423", downregulated = "#2878B5")) +
  labs(title = "Differential Expression Summary", x = "Comparison", y = "Number of Significant Genes", fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plots_dir, "deg_summary_barplot.pdf"), deg_barplot, width = 10, height = 6)

total_deg_plot <- ggplot(comparison_summary_df, aes(reorder(comparison, total_significant), total_significant)) +
  geom_col(fill = "grey40") + coord_flip() +
  labs(title = "Total Significant Genes by Comparison", x = "", y = "Significant Genes") + theme_classic()
ggsave(file.path(plots_dir, "total_significant_genes.pdf"), total_deg_plot, width = 8, height = 6)

result_files <- list.files(results_dir, pattern = "\\.csv$", full.names = TRUE)
combined_deg_table <- map_dfr(result_files, function(file) {
  dat <- read.csv(file, stringsAsFactors = FALSE)
  dat$comparison <- sub("\\.csv$", "", basename(file))
  dat
})
write.csv(combined_deg_table, file.path(output_dir, "all_comparisons_combined.csv"), row.names = FALSE)

vst_matrix <- assay(vsd)
gene_variance <- apply(vst_matrix, 1, var)
variance_table <- tibble(gene_id = names(gene_variance), variance = as.numeric(gene_variance)) |>
  arrange(desc(variance))
write.csv(variance_table, file.path(output_dir, "gene_variance_ranking.csv"), row.names = FALSE)
top_variable_genes <- slice_head(variance_table, n = min(100, nrow(variance_table)))
write.csv(top_variable_genes, file.path(output_dir, "top100_variable_genes.csv"), row.names = FALSE)
heatmap_matrix <- vst_matrix[top_variable_genes$gene_id, , drop = FALSE]
pdf(file.path(qc_dir, "qc_top100_variable_genes_heatmap.pdf"), width = 10, height = 10)
pheatmap(heatmap_matrix, scale = "row", annotation_col = annotation_df, show_rownames = FALSE, main = "Top 100 Most Variable Genes")
dev.off()

analysis_timestamp <- Sys.time()
write.csv(
  tibble(
    metric = c("analysis_timestamp", "samples", "genes_tested", "comparisons", "alpha", "log2fc_threshold"),
    value = c(as.character(analysis_timestamp), ncol(dds), nrow(dds), nrow(comparison_summary_df), alpha, log2fc_threshold)
  ),
  file.path(output_dir, "analysis_report.csv"), row.names = FALSE
)
write.csv(
  tibble(
    package = required_packages,
    version = vapply(required_packages, function(x) as.character(packageVersion(x)), character(1))
  ),
  file.path(output_dir, "software_versions.csv"), row.names = FALSE
)
write.csv(
  metadata |> mutate(size_factor = sizeFactors(dds)) |> arrange(genotype, replicate),
  file.path(output_dir, "sample_summary.csv"), row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

expected_outputs <- c(
  file.path(output_dir, "raw_counts.csv"),
  file.path(output_dir, "normalized_counts.csv"),
  file.path(output_dir, "vst_counts.csv"),
  file.path(output_dir, "comparison_summary.csv"),
  file.path(qc_dir, "qc_pca.pdf")
)
missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
if (length(missing_outputs)) stop("Missing outputs: ", paste(missing_outputs, collapse = ", "), call. = FALSE)

writeLines(
  c("DESeq2 analysis completed successfully.", paste("Completed:", analysis_timestamp)),
  file.path(output_dir, "analysis_complete.flag")
)
write.csv(
  tibble(completed = TRUE, completion_time = as.character(Sys.time()), samples = ncol(dds), genes_tested = nrow(dds), comparisons = nrow(comparison_summary_df)),
  file.path(output_dir, "run_completion.csv"), row.names = FALSE
)

cat("\nRNA-SEQ ANALYSIS COMPLETE\n")
cat("Samples:", ncol(dds), "\n")
cat("Genes tested:", nrow(dds), "\n")
cat("Comparisons:", nrow(comparison_summary_df), "\n")
cat("Output directory:", output_dir, "\n")
quit(save = "no", status = 0)
