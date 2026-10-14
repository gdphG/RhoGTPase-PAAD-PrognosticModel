###############################################################################
# Step 5: LASSO-Cox Regression (CORE MODEL BUILDING)
# - 10-fold cross-validation for lambda selection
# - Bootstrap resampling (200 iterations) for Optimism-corrected AUC
# - Permutation test (500 iterations) for statistical significance
###############################################################################

library(survival)
library(glmnet)
library(timeROC)

# --- Load data ---
risk_data <- read.table("data/risk_scorecd.txt", header = TRUE, sep = "\t")
exp_data <- read.table("data/timeexpcd.txt", header = TRUE, sep = "\t", 
                        row.names = 1, check.names = FALSE)

# Merge: expression + risk scores
row.names(risk_data) <- risk_data[, 1]
merged <- merge(exp_data, risk_data, by = "row.names")
row.names(merged) <- merged[, 1]
merged <- merged[, -1]

# --- Prepare LASSO input ---
# Predictor matrix (gene expression)
gene_cols <- setdiff(colnames(merged), c("futime", "fustat", "riskScore", "Risk"))
v1 <- t(as.matrix(merged[, gene_cols]))  # genes × samples
v1 <- apply(v1, 2, as.numeric)

# Response (survival)
v2 <- Surv(merged$futime, merged$fustat)

# --- LASSO-Cox with cross-validation ---
set.seed(42)
myfit <- glmnet(v1, v2, family = "cox", alpha = 1)  # alpha=1 for LASSO
myfit2 <- cv.glmnet(v1, v2, family = "cox", nfolds = 10)

# Extract non-zero coefficient genes
lambda_min <- myfit2$lambda.min
coef_min <- coef(myfit, s = lambda_min)
selected_genes <- rownames(coef_min)[which(coef_min != 0)]

cat("LASSO selected", length(selected_genes), "genes:\n")
cat(paste(selected_genes, collapse = ", "), "\n")

# --- Plot LASSO results ---
pdf("figures/lasso_lambda_path.pdf", width = 8, height = 6)
plot(myfit2, xlab = "log(lambda)", ylab = "Cross-validation deviance")
abline(v = log(lambda_min), col = "red", lty = 2)
dev.off()

# --- Save LASSO-selected genes ---
lasso_result <- data.frame(
  gene = selected_genes,
  coefficient = as.numeric(coef_min[which(coef_min != 0)])
)
write.table(lasso_result, "data/lasso_selected_genes.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

###############################################################################
# Bootstrap Analysis: Optimism-corrected AUC (200 iterations)
###############################################################################

n_boot <- 200
auc_original <- numeric(n_boot)
auc_optimism <- numeric(n_boot)

set.seed(123)
for (i in 1:n_boot) {
  # Resample with replacement
  boot_idx <- sample(1:nrow(merged), nrow(merged), replace = TRUE)
  boot_data <- merged[boot_idx, ]
  
  # Fit LASSO on bootstrap sample
  v1_boot <- t(as.matrix(boot_data[, selected_genes]))
  v1_boot <- apply(v1_boot, 2, as.numeric)
  v2_boot <- Surv(boot_data$futime, boot_data$fustat)
  
  fit_boot <- cv.glmnet(v1_boot, v2_boot, family = "cox", nfolds = 10)
  coef_boot <- coef(fit_boot, s = "lambda.min")
  
  # Calculate risk score on bootstrap sample
  risk_boot <- as.numeric(v1_boot %*% coef_boot[-1])  # exclude intercept
  
  # Time-dependent ROC (1-year)
  roc_boot <- timeROC(T = boot_data$futime, delta = boot_data$fustat,
                       marker = risk_boot, cause = 1, weighting = "marginal",
                       times = 365)
  auc_original[i] <- roc_boot$AUC[2]
  
  # Apply to original data for optimism estimation
  v1_orig <- t(as.matrix(merged[, selected_genes]))
  v1_orig <- apply(v1_orig, 2, as.numeric)
  risk_orig <- as.numeric(v1_orig %*% coef_boot[-1])
  
  roc_orig <- timeROC(T = merged$futime, delta = merged$fustat,
                       marker = risk_orig, cause = 1, weighting = "marginal",
                       times = 365)
  auc_optimism[i] <- roc_orig$AUC[2]
}

# Optimism-corrected AUC
optimism <- mean(auc_optimism) - mean(auc_original)
auc_corrected <- mean(auc_original) - optimism
cat("Optimism-corrected AUC (1-year):", round(auc_corrected, 4), "\n")

###############################################################################
# Permutation Test (500 iterations)
###############################################################################

n_perm <- 500
auc_perm <- numeric(n_perm)

set.seed(456)
for (i in 1:n_perm) {
  # Permute survival outcomes
  perm_surv_time <- sample(merged$futime)
  perm_surv_status <- merged$fustat
  
  v1_perm <- t(as.matrix(merged[, selected_genes]))
  v1_perm <- apply(v1_perm, 2, as.numeric)
  v2_perm <- Surv(perm_surv_time, perm_surv_status)
  
  fit_perm <- cv.glmnet(v1_perm, v2_perm, family = "cox", nfolds = 10)
  coef_perm <- coef(fit_perm, s = "lambda.min")
  risk_perm <- as.numeric(v1_perm %*% coef_perm[-1])
  
  roc_perm <- timeROC(T = perm_surv_time, delta = perm_surv_status,
                       marker = risk_perm, cause = 1, weighting = "marginal",
                       times = 365)
  auc_perm[i] <- roc_perm$AUC[2]
}

# Permutation P-value
p_perm <- mean(auc_perm >= mean(auc_original))
cat("Permutation P-value:", p_perm, "\n")

# --- Save results ---
bootstrap_result <- data.frame(
  metric = c("Mean_AUC_original", "Mean_AUC_optimism", "Optimism", 
             "Corrected_AUC", "Permutation_P"),
  value = c(mean(auc_original), mean(auc_optimism), optimism, 
            auc_corrected, p_perm)
)
write.table(bootstrap_result, "data/lasso_validation_results.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
