###############################################################################
# Step 10: Clinical Feature ROC Comparison
# Compare risk score AUC vs. individual clinical features at 1/2/3/5 years
###############################################################################

library(timeROC)

# --- Load data ---
rt <- read.table("data/riskcli_subset.txt", header = TRUE, sep = "\t")
colnames(rt)[1:3] <- c("futime", "fustat", "riskScore")

# --- Function to plot clinical ROC at a given time point ---
plot_clinical_ROC <- function(rt, predictTime, output_file) {
  risk <- rt[, 1:3]
  bioCol <- rainbow(ncol(rt) - 1, s = 0.8, v = 0.8)
  aucText <- c()
  
  pdf(file = output_file, width = 5.5, height = 5.5)
  
  # Risk score ROC
  ROC_rt <- timeROC(T = risk$futime, delta = risk$fustat,
                     marker = risk$riskScore, cause = 1,
                     weighting = "aalen", times = predictTime, ROC = TRUE)
  plot(ROC_rt, time = predictTime, col = bioCol[1], title = FALSE, lwd = 2)
  aucText <- c(paste0("Riskscore, AUC=", sprintf("%.3f", ROC_rt$AUC[2])))
  abline(0, 1)
  
  # Clinical features ROC
  for (i in 4:ncol(rt)) {
    ROC_rt <- timeROC(T = rt$futime, delta = rt$fustat,
                       marker = rt[, i], cause = 1,
                       weighting = "aalen", times = predictTime, ROC = TRUE)
    plot(ROC_rt, time = predictTime, col = bioCol[i - 2], 
         title = FALSE, lwd = 2, add = TRUE)
    aucText <- c(aucText, paste0(colnames(rt)[i], ", AUC=", 
                                  sprintf("%.3f", ROC_rt$AUC[2])))
  }
  
  legend("bottomright", aucText, lwd = 2, bty = "n", 
         col = bioCol[1:(ncol(rt) - 1)], cex = 0.7)
  dev.off()
}

# --- Generate ROC plots for each time point ---
plot_clinical_ROC(rt, 365, "figures/cliROC_1year.pdf")
plot_clinical_ROC(rt, 730, "figures/cliROC_2year.pdf")
plot_clinical_ROC(rt, 1095, "figures/cliROC_3year.pdf")
plot_clinical_ROC(rt, 1825, "figures/cliROC_5year.pdf")
