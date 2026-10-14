###############################################################################
# Step 16: Drug Sensitivity Analysis
# Using pRRophetic algorithm with GDSC database
# FDR correction for multiple testing
###############################################################################

library(pRRophetic)
library(dplyr)

# --- Load data ---
geneexp <- read.table("data/mRNA1.txt", header = TRUE, sep = "\t", 
                       check.names = FALSE)
risk_data <- read.table("data/risk_scores_final.txt", header = TRUE, sep = "\t")

# Prepare expression matrix
expr_matrix <- as.matrix(geneexp[, -1])
rownames(expr_matrix) <- geneexp$id
colnames(expr_matrix) <- substr(colnames(expr_matrix), 1, 12)

# Assign risk groups
risk_data$risk_level <- ifelse(risk_data$risk_score > median(risk_data$risk_score), 
                                "High", "Low")

# --- Drug sensitivity prediction ---
# Get common samples between expression and risk data
high_samples <- risk_data$id[risk_data$risk_level == "High"]
low_samples <- risk_data$id[risk_data$risk_level == "Low"]

# pRRophetic analysis for each drug
drugs <- c("Gemcitabine", "Paclitaxel", "Cisplatin", "5-Fluorouracil",
           "Erlotinib", "Sorafenib", "Olaparib", "Nab-paclitaxel")

results <- data.frame()

for (drug in drugs) {
  tryCatch({
    # Predict IC50 for high-risk samples
    high_expr <- expr_matrix[, high_samples, drop = FALSE]
    pRR.high <- pRRophetic.predict(inputGeneExp = high_expr, 
                                    drugName = drug,
                                    tissueType = "all")
    
    # Predict IC50 for low-risk samples
    low_expr <- expr_matrix[, low_samples, drop = FALSE]
    pRR.low <- pRRophetic.predict(inputGeneExp = low_expr,
                                   drugName = drug,
                                   tissueType = "all")
    
    # Wilcoxon test
    test <- wilcox.test(pRR.high, pRR.low)
    
    # Effect size (Cohen's d)
    pooled_sd <- sqrt((var(pRR.high) + var(pRR.low)) / 2)
    cohens_d <- (mean(pRR.high) - mean(pRR.low)) / pooled_sd
    
    results <- rbind(results, data.frame(
      Drug = drug,
      Mean_IC50_High = mean(pRR.high),
      Mean_IC50_Low = mean(pRR.low),
      P_value = test$p.value,
      Cohens_d = cohens_d
    ))
  }, error = function(e) {
    cat("Error for", drug, ":", conditionMessage(e), "\n")
  })
}

# --- FDR correction ---
results$FDR <- p.adjust(results$P_value, method = "BH")
write.table(results, "data/drug_sensitivity_results.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Visualization ---
results$Significant <- ifelse(results$FDR < 0.05, "Yes", "No")

pdf("figures/drug_sensitivity_volcano.pdf", width = 8, height = 6)
ggplot(results, aes(x = Cohens_d, y = -log10(FDR), 
                     color = Significant, label = Drug)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  xlab("Cohen's d (High - Low)") +
  ylab("-log10(FDR)") +
  ggtitle("Drug Sensitivity: High vs. Low Risk Group")
dev.off()
