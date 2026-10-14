###############################################################################
# Step 13: CIBERSORT Immune Infiltration + Risk Score Correlation
# Correlate 22 immune cell types with risk score
###############################################################################

library(ggpubr)
library(dplyr)
library(readr)

# --- Load CIBERSORT results ---
cibersort <- read_delim("data/CIBERSORT-Results.txt", "\t", 
                         escape_double = FALSE, trim_ws = TRUE)
risk_data <- read.table("data/risk_scores_final.txt", header = TRUE, sep = "\t")
risk_data <- risk_data[, -c(2:4)]
colnames(risk_data)[2] <- "Riskscore"

# --- Merge ---
cibersort$id <- substr(cibersort$CellType, 1, 12)
merged <- merge(cibersort, risk_data, by = "id")

# --- Correlation analysis ---
# Extract cell type proportions (columns 2:23)
cell_types <- colnames(merged)[2:23]

results <- data.frame()
for (ct in cell_types) {
  cor_test <- cor.test(merged[, ct], merged$Riskscore, method = "spearman")
  results <- rbind(results, data.frame(
    CellType = ct,
    rho = cor_test$estimate,
    P_value = cor_test$p.value
  ))
}
results$FDR <- p.adjust(results$P_value, method = "BH")
write.table(results, "data/immune_correlation_results.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Visualization: Box plots by risk group ---
sdata <- stack(merged[, cell_types])
sdata$Riskscore <- ifelse(merged$Riskscore > median(merged$Riskscore), "High", "Low")

pdf("figures/immune_composition_by_riskgroup.pdf", width = 8, height = 6)
ggboxplot(sdata, x = "ind", y = "values",
          color = "black", fill = "Riskscore",
          xlab = "", ylab = "Cell composition",
          main = "TME Cell composition by Risk Score") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()
