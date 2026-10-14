###############################################################################
# Step 1: Differential Gene Expression Analysis
# Wilcoxon rank-sum test with FDR correction (Benjamini-Hochberg)
# Criteria: |log2FC| > 1 & FDR < 0.05
###############################################################################

library(limma)

# --- Load data ---
# geneexp: expression matrix (genes × samples), first column = gene symbols
# Input format: tab-delimited, row 1 = header (sample IDs), col 1 = gene symbols
geneexp <- read.table("data/mRNA1.txt", header = TRUE, sep = "\t", 
                       check.names = FALSE, row.names = 1)

# --- Define groups ---
# Extract TCGA barcode: 14th-15th character indicates sample type
# 01 = tumor, 11 = normal
group <- sapply(strsplit(colnames(geneexp), "-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
tumor_idx <- which(group == 1)   # Tumor samples
normal_idx <- which(group == 2)  # Normal samples (recoded from 11->2)

group1 <- length(tumor_idx)
group2 <- length(normal_idx)

# Prepare combined data for Wilcoxon test
newgene <- geneexp

# --- Wilcoxon rank-sum test for each gene ---
results <- data.frame()

for (a in row.names(newgene)) {
  gene_expr <- as.numeric(newgene[a, ])
  mydata <- data.frame(
    geneexp = c(gene_expr[tumor_idx], gene_expr[normal_idx]),
    mygroup = c(rep(1, group1), rep(2, group2))
  )
  
  # Calculate log2FC
  group1mean <- mean(gene_expr[tumor_idx])
  group2mean <- mean(gene_expr[normal_idx])
  logFC <- log2(group2mean) - log2(group1mean)
  
  # Wilcoxon test
  wilcoxTest <- wilcox.test(geneexp ~ mygroup, data = mydata)
  
  results <- rbind(results, data.frame(
    gene = a,
    logFC = logFC,
    P_value = wilcoxTest$p.value
  ))
}

# --- FDR correction (Benjamini-Hochberg) ---
results$FDR <- p.adjust(results$P_value, method = "BH")

# --- Filter significant DEGs ---
deg <- results[abs(results$logFC) > 1 & results$FDR < 0.05, ]
write.table(deg, "data/differential_genes.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)

cat("Number of differentially expressed genes:", nrow(deg), "\n")
