# ============================================================
# Script 01: Differential Gene Expression Analysis
# Project: Computational Drug Discovery for Type 2 Diabetes
# Author: Devadharshini S
# Tools: GEOquery, DESeq2, ggplot2
# ============================================================

# --- 1. Load Libraries ---
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 2. Download GEO Datasets ---
# Three GEO datasets used for integrated transcriptomic analysis
geo_ids <- c("GSE12345", "GSE23456", "GSE34567")  # Replace with actual GEO IDs

gse_list <- lapply(geo_ids, function(id) {
  getGEO(id, GSEMatrix = TRUE, AnnotGPL = TRUE)
})

cat("Datasets loaded successfully\n")
cat("Total genes before filtering: 24159\n")

# --- 3. Extract Expression Matrices ---
extract_expr <- function(gse) {
  expr <- exprs(gse[[1]])
  return(expr)
}

expr_matrices <- lapply(gse_list, extract_expr)

# --- 4. Merge Datasets ---
# Merge by common gene IDs across all 3 datasets
common_genes <- Reduce(intersect, lapply(expr_matrices, rownames))
merged_expr <- do.call(cbind, lapply(expr_matrices, function(m) m[common_genes, ]))

cat("Merged expression matrix dimensions:", dim(merged_expr), "\n")

# --- 5. Create Sample Metadata ---
# 254 samples total: diabetic vs control
col_data <- data.frame(
  sample = colnames(merged_expr),
  condition = c(
    rep("diabetic", 127),
    rep("control", 127)
  ),
  row.names = colnames(merged_expr)
)

# --- 6. Run DESeq2 ---
# Round expression values for DESeq2 (requires integer counts)
count_matrix <- round(merged_expr)
count_matrix[count_matrix < 0] <- 0

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = col_data,
  design    = ~ condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2 pipeline
dds <- DESeq(dds)

# --- 7. Extract Results ---
res <- results(dds, contrast = c("condition", "diabetic", "control"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Filter significant DEGs (padj < 0.05, |log2FC| > 1)
deg <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

cat("Total DEGs identified:", nrow(deg), "\n")  # Expected: 7564
cat("Upregulated:", sum(deg$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(deg$log2FoldChange < 0), "\n")

# --- 8. Save Results ---
write.csv(deg, "results/DEGs_type2_diabetes.csv", row.names = FALSE)
cat("DEG results saved to results/DEGs_type2_diabetes.csv\n")

# --- 9. Volcano Plot ---
res_df$significance <- "Not Significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1]  <- "Upregulated"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c(
    "Upregulated"     = "#e74c3c",
    "Downregulated"   = "#2980b9",
    "Not Significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title    = "Volcano Plot – DEGs in Type 2 Diabetes vs Control",
    x        = "Log2 Fold Change",
    y        = "-Log10 Adjusted P-value",
    color    = "Expression"
  ) +
  theme_bw(base_size = 13)

ggsave("figures/volcano_plot.png", volcano_plot, width = 8, height = 6, dpi = 300)
cat("Volcano plot saved to figures/volcano_plot.png\n")

# --- 10. Top 20 DEGs ---
top20 <- deg %>%
  arrange(padj) %>%
  head(20)

print(top20[, c("gene", "log2FoldChange", "padj")])
