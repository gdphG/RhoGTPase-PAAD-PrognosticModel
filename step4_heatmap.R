###############################################################################
# Step 4: Heatmap Visualization of Significant Genes
# Using pheatmap package
###############################################################################

library(pheatmap)

# --- Load expression data and gene list ---
geneexp <- read.table("data/mRNA1.txt", header = TRUE, sep = "\t", 
                       check.names = FALSE, row.names = 1)

overlap_genes <- read.table("data/overlap_genes.txt", header = FALSE)[, 1]

# Filter expression matrix
exp_subset <- geneexp[rownames(geneexp) %in% overlap_genes, ]

# Log2 transform (if not already)
exp_log <- log2(exp_subset + 1)

# Scale (z-score)
exp_scaled <- t(scale(t(exp_log)))

# --- Annotation ---
# Tumor vs Normal annotation
group <- sapply(strsplit(colnames(exp_scaled), "-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
group <- gsub("2", "1", group)
annotation_col <- data.frame(
  Type = ifelse(group == 1, "Tumor", "Normal"),
  row.names = colnames(exp_scaled)
)

# --- Plot ---
pdf("figures/heatmap_significant_genes.pdf", width = 12, height = 10)
pheatmap(exp_scaled,
         annotation_col = annotation_col,
         show_colnames = FALSE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Expression of Significant Rho GTPase Genes")
dev.off()
