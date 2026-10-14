###############################################################################
# Step 19: External Validation
# Validate 22-gene model in independent cohorts:
# GSE28735, GSE62452, ICGC-PACA-CA, GSE57495, GSE78220, IMvigor210
###############################################################################

library(survival)
library(survminer)
library(timeROC)

# --- 22-gene signature with coefficients ---
signature_genes <- c("ARAP1", "ARHGAP18", "ARHGAP23", "ARHGAP29", "ARHGAP44",
                     "ARHGEF11", "ARHGEF26", "CALM1", "CENPT", "CTNNA1",
                     "CTTN", "DEPDC1B", "DVL1", "GDI1", "INPP5B", "IQGAP1",
                     "NCKAP1", "NUP37", "PAK2", "RHOBTB2", "RHOV", "SPDL1")

# Load coefficients from multivariate Cox (Step 6)
multi_cox <- read.table("data/multivariate_cox_results.txt", header = TRUE, sep = "\t")
coef_vec <- setNames(multi_cox$coef, multi_cox$gene)

# Median risk score from training cohort (TCGA-PAAD)
training_median <- 0.0  # Replace with actual value from step6 output

# --- Function for validation ---
validate_cohort <- function(exp_file, clinical_file, cohort_name, 
                             output_dir = "figures/validation/") {
  
  # Load data
  exp <- read.table(exp_file, header = TRUE, sep = "\t", 
                     check.names = FALSE, row.names = 1)
  cli <- read.table(clinical_file, header = TRUE, sep = "\t")
  
  # Check gene availability
  available_genes <- intersect(signature_genes, rownames(exp))
  cat(cohort_name, ": Available genes =", length(available_genes), "/", 
      length(signature_genes), "\n")
  
  # Calculate risk score
  common_samples <- intersect(colnames(exp), cli$id)
  exp_sub <- exp[, common_samples]
  cli <- cli[match(common_samples, cli$id), ]
  
  # Risk score = sum(coef * expression)
  risk_score <- numeric(length(common_samples))
  for (gene in available_genes) {
    if (gene %in% names(coef_vec)) {
      risk_score <- risk_score + coef_vec[gene] * as.numeric(exp_sub[gene, ])
    }
  }
  
  risk_level <- ifelse(risk_score > training_median, "High", "Low")
  
  # Survival analysis
  df <- data.frame(
    survival_time = cli$survival_time,
    status = cli$status,
    risk_score = risk_score,
    risk_level = risk_level
  )
  
  fit <- survfit(Surv(survival_time, status) ~ risk_level, data = df)
  pval <- survdiff(Surv(survival_time, status) ~ risk_level, data = df)$pvalue
  
  # KM plot
  pdf(paste0(output_dir, "validation_KM_", cohort_name, ".pdf"), 
      width = 8, height = 7)
  print(ggsurvplot(fit, data = df, 
                    conf.int = TRUE, risk.table = TRUE,
                    palette = c("#E41A1C", "#377EB8"),
                    title = paste0(cohort_name, " (P = ", 
                                   format(pval, digits = 3), ")"),
                    legend.labs = c("High", "Low")))
  dev.off()
  
  # Time-dependent ROC
  tryCatch({
    roc <- timeROC(T = df$survival_time, delta = df$status,
                    marker = df$risk_score, cause = 1,
                    weighting = "marginal", times = c(365, 1095, 1825))
    
    pdf(paste0(output_dir, "validation_ROC_", cohort_name, ".pdf"),
        width = 10, height = 3.5)
    par(mfrow = c(1, 3))
    for (t in c(365, 1095, 1825)) {
      label <- paste0(t/365, "-year")
      plot(roc, time = t, col = "red", lwd = 2,
           title = paste0(label, " AUC=", sprintf("%.3f", roc$AUC[which(roc$times == t)])))
      abline(0, 1, col = "gray", lty = 2)
    }
    dev.off()
  }, error = function(e) {
    cat("ROC error for", cohort_name, ":", conditionMessage(e), "\n")
  })
  
  return(data.frame(
    Cohort = cohort_name,
    N = nrow(df),
    N_High = sum(risk_level == "High"),
    N_Low = sum(risk_level == "Low"),
    P_value = pval
  ))
}

# --- Run validation for each cohort ---
dir.create("figures/validation", showWarnings = FALSE)

results_list <- list()

# GSE28735
results_list[[1]] <- validate_cohort("data/GSE28735/expression.txt",
                                      "data/GSE28735/clinical.txt", "GSE28735")

# GSE62452
results_list[[2]] <- validate_cohort("data/GSE62452/expression.txt",
                                      "data/GSE62452/clinical.txt", "GSE62452")

# ICGC-PACA-CA
results_list[[3]] <- validate_cohort("data/ICGC/expression.txt",
                                      "data/ICGC/clinical.txt", "ICGC-PACA-CA")

# GSE57495 (biopsy)
results_list[[4]] <- validate_cohort("data/GSE57495/expression.txt",
                                      "data/GSE57495/clinical.txt", "GSE57495")

# GSE78220 (immunotherapy)
results_list[[5]] <- validate_cohort("data/GSE78220/expression.txt",
                                      "data/GSE78220/clinical.txt", "GSE78220")

# IMvigor210 (immunotherapy)
results_list[[6]] <- validate_cohort("data/IMvigor210/expression.txt",
                                      "data/IMvigor210/clinical.txt", "IMvigor210")

# --- Summary table ---
validation_summary <- do.call(rbind, results_list)
write.table(validation_summary, "data/validation_summary.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Validation Summary ===\n")
print(validation_summary)
