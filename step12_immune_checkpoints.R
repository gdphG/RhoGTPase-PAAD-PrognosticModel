###############################################################################
# Step 12: Immune Checkpoint Gene Expression Analysis
# Compare expression of 8 immune checkpoint genes between High/Low risk
###############################################################################

library(ggpubr)
library(dplyr)

# --- Immune checkpoint genes ---
icp_genes <- c("CD274", "CTLA4", "HAVCR2", "LAG3", 
               "PDCD1", "PDCD1LG2", "SIGLEC15", "TIGIT")

# --- Load data ---
geneexp <- read.table("data/mRNA1.txt", header = TRUE, sep = "\t", 
                       check.names = FALSE)
risk_data <- read.table("data/risk_scores_final.txt", header = TRUE, sep = "\t")

# Extract ICP gene expression
icp_expr <- geneexp[geneexp$id %in% icp_genes, ]
icp_expr <- icp_expr[, -c(2:5)]  # Remove non-sample columns
colnames(icp_expr) <- substr(colnames(icp_expr), 1, 12)  # Truncate to 12 chars

# Map sample IDs to risk groups
id_to_risk <- setNames(risk_data$risk_level, risk_data$id)
colnames(icp_expr) <- id_to_risk[colnames(icp_expr)]

# Transpose for analysis (samples × genes)
icp_t <- t(icp_expr)
icp_df <- as.data.frame(icp_t)
icp_df$RiskGroup <- rownames(icp_df)

# --- Statistical comparison (Wilcoxon test) ---
results <- data.frame()
for (gene in icp_genes) {
  high_vals <- icp_df[icp_df$RiskGroup == "High", gene]
  low_vals <- icp_df[icp_df$RiskGroup == "Low", gene]
  
  test <- wilcox.test(high_vals, low_vals)
  results <- rbind(results, data.frame(
    gene = gene,
    mean_High = mean(high_vals),
    mean_Low = mean(low_vals),
    P_value = test$p.value
  ))
}

# FDR correction
results$FDR <- p.adjust(results$P_value, method = "BH")
write.table(results, "data/immune_checkpoint_comparison.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Visualization ---
# Prepare long-format data for plotting
icp_long <- icp_df %>%
  tidyr::pivot_longer(cols = icp_genes, names_to = "Gene", values_to = "Expression")

pdf("figures/immune_checkpoint_boxplot.pdf", width = 10, height = 8)
ggboxplot(icp_long, x = "Gene", y = "Expression", 
          fill = "RiskGroup", palette = c("#E41A1C", "#377EB8"),
          add = "jitter") +
  stat_compare_means(aes(group = RiskGroup), label = "p.signif", 
                     method = "wilcox.test") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Immune Checkpoint Gene Expression: High vs. Low Risk")
dev.off()
