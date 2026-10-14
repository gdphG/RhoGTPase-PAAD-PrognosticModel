###############################################################################
# Step 17: Immunotherapy Response Prediction
# IPS scoring + TIDE algorithm for immunotherapy response
###############################################################################

# --- IPS (Immunophenoscore) Analysis ---
# Download from tcia.at (The Cancer Immunome Atlas)
# IPS data: http://tcia.at/immunophenogram/

library(dplyr)

# IPS values (pre-downloaded)
ips_data <- read.table("data/ips_PAAD.txt", header = TRUE, sep = "\t", 
                        check.names = FALSE)
risk_data <- read.table("data/risk_scores_final.txt", header = TRUE, sep = "\t")

# Merge
ips_data$id <- substr(ips_data$id, 1, 12)
risk_data$id_short <- substr(risk_data$id, 1, 12)
merged <- merge(ips_data, risk_data, by.x = "id", by.y = "id_short")

# IPS categories: ACT, EC, MHC, SCK
for (category in c("ACT", "EC", "MHC", "SCK")) {
  if (category %in% colnames(merged)) {
    high_vals <- merged[merged$risk_level == "High", category]
    low_vals <- merged[merged$risk_level == "Low", category]
    test <- wilcox.test(high_vals, low_vals)
    cat(category, ": P =", test$p.value, "\n")
  }
}

# --- TIDE Algorithm ---
# Download from: http://tide.dfci.harvard.edu/
# Upload expression matrix and get TIDE prediction

tide_results <- read.table("data/tide_PAAD_results.txt", header = TRUE, sep = "\t")

# Merge with risk data
tide_results$id <- substr(tide_results$id, 1, 12)
merged_tide <- merge(tide_results, risk_data[, c("id", "risk_level")], by = "id")

# Compare TIDE scores
high_tide <- merged_tide[merged_tide$risk_level == "High", "TIDE"]
low_tide <- merged_tide[merged_tide$risk_level == "Low", "TIDE"]
wilcox.test(high_tide, low_tide)

# --- Response rate comparison ---
# TIDE responder: non-responder = TRUE
response_table <- table(
  Response = merged_tide$Responder,
  RiskGroup = merged_tide$risk_level
)
chisq.test(response_table)

write.table(response_table, "data/immunotherapy_response_table.txt",
            sep = "\t", quote = TRUE)

# --- Visualization ---
pdf("figures/TIDE_score_comparison.pdf", width = 5, height = 5)
boxplot(TIDE ~ risk_level, data = merged_tide,
        col = c("#377EB8", "#E41A1C"),
        main = "TIDE Score: High vs. Low Risk",
        ylab = "TIDE score")
stripchart(TIDE ~ risk_level, data = merged_tide, 
           vertical = TRUE, method = "jitter", 
           add = TRUE, pch = 1, col = "gray")
dev.off()
