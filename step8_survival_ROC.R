###############################################################################
# Step 8: Kaplan-Meier Survival Curves + Time-dependent ROC
# 1-, 2-, 3-, 5-year time-dependent ROC analysis
###############################################################################

library(survival)
library(survminer)
library(timeROC)

# --- Load data ---
rt <- read.table("data/risk_scores_final.txt", header = TRUE, sep = "\t")
clinical <- read.table("data/clinical.txt", header = TRUE, sep = "\t")
clinical <- clinical[complete.cases(clinical[, 2]), ]

# Merge
data1 <- merge(rt, clinical[, c("id", "survival_time", "status")], by.x = "id", by.y = "id")

# --- Kaplan-Meier Survival Analysis ---
fit <- survfit(Surv(survival_time, status) ~ risk_level, data = data1)

# KM plot with risk table
surPlot <- ggsurvplot(
  fit,
  data = data1,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#E41A1C", "#377EB8"),  # High=red, Low=blue
  xlab = "Time (days)",
  ylab = "Survival probability",
  title = "Overall Survival: High vs. Low Risk",
  legend.title = "Risk group",
  legend.labs = c("High", "Low"),
  pval = TRUE,
  pval.size = 5,
  ggtheme = theme_bw()
)

pdf("figures/km_survival.pdf", width = 8, height = 7)
print(surPlot)
dev.off()

# Log-rank test
survdiff(Surv(survival_time, status) ~ risk_level, data = data1)

# --- Time-dependent ROC Analysis ---
TCGA <- data1

# 1-year ROC
ROC_1y <- timeROC(T = TCGA$survival_time, delta = TCGA$status,
                   marker = TCGA$risk_score, cause = 1,
                   weighting = "marginal", times = 365)

# 2-year ROC
ROC_2y <- timeROC(T = TCGA$survival_time, delta = TCGA$status,
                   marker = TCGA$risk_score, cause = 1,
                   weighting = "marginal", times = 1095)

# 3-year ROC
ROC_3y <- timeROC(T = TCGA$survival_time, delta = TCGA$status,
                   marker = TCGA$risk_score, cause = 1,
                   weighting = "marginal", times = 1825)

# 5-year ROC
ROC_5y <- timeROC(T = TCGA$survival_time, delta = TCGA$status,
                   marker = TCGA$risk_score, cause = 1,
                   weighting = "marginal", times = 1825)

# --- Plot ROC curves ---
pdf("figures/time_dependent_ROC.pdf", width = 10, height = 10)
par(mfrow = c(2, 2))

# 1-year
plot(ROC_1y, time = 365, col = "red", lwd = 2, 
     title = paste0("1-year ROC (AUC = ", sprintf("%.3f", ROC_1y$AUC[2]), ")"))
abline(0, 1, col = "gray", lty = 2)

# 2-year
plot(ROC_2y, time = 1095, col = "blue", lwd = 2,
     title = paste0("2-year ROC (AUC = ", sprintf("%.3f", ROC_2y$AUC[2]), ")"))
abline(0, 1, col = "gray", lty = 2)

# 3-year
plot(ROC_3y, time = 1825, col = "green", lwd = 2,
     title = paste0("3-year ROC (AUC = ", sprintf("%.3f", ROC_3y$AUC[2]), ")"))
abline(0, 1, col = "gray", lty = 2)

# 5-year
plot(ROC_5y, time = 1825, col = "purple", lwd = 2,
     title = paste0("5-year ROC (AUC = ", sprintf("%.3f", ROC_5y$AUC[2]), ")"))
abline(0, 1, col = "gray", lty = 2)

dev.off()

# --- Save ROC results ---
roc_summary <- data.frame(
  Time = c("1-year", "2-year", "3-year", "5-year"),
  AUC = c(ROC_1y$AUC[2], ROC_2y$AUC[2], ROC_3y$AUC[2], ROC_5y$AUC[2])
)
write.table(roc_summary, "data/roc_auc_summary.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
