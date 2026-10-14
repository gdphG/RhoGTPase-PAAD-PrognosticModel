###############################################################################
# Step 3: Venn Diagram - Gene Intersection
# Identify overlapping genes between differential expression and 
# univariate Cox significant genes
###############################################################################

library(VennDiagram)

# --- Load gene lists ---
# DEGs from Step 1
deg_genes <- read.table("data/differential_genes.txt", header = TRUE, sep = "\t")$gene

# Significant Cox genes from Step 2
cox_genes <- read.table("data/univariate_cox_significant.txt", header = TRUE, sep = "\t")$gene

# --- Venn diagram ---
pdf("figures/venn_diagram.pdf", width = 6, height = 6)
venn.plot <- venn.diagram(
  x = list(
    "Differential\nExpression" = deg_genes,
    "Cox Regression" = cox_genes
  ),
  filename = NULL,
  col = "transparent",
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.fontface = "bold",
  lwd = 2,
  main = "Intersection of DEGs and Cox-significant genes"
)
grid.draw(venn.plot)
dev.off()

# Extract intersection genes
overlap_genes <- intersect(deg_genes, cox_genes)
write.table(overlap_genes, "data/overlap_genes.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Overlapping genes:", length(overlap_genes), "\n")
