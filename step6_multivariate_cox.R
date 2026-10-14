###############################################################################
# Step 6: Multivariate Cox Regression + Risk Score + Forest Plot
# Constructs the final 22-gene prognostic model
###############################################################################

library(survival)
library(forestplot)

# --- Define 22 signature genes ---
signature_genes <- c("ARAP1", "ARHGAP18", "ARHGAP23", "ARHGAP29", "ARHGAP44",
                     "ARHGEF11", "ARHGEF26", "CALM1", "CENPT", "CTNNA1",
                     "CTTN", "DEPDC1B", "DVL1", "GDI1", "INPP5B", "IQGAP1",
                     "NCKAP1", "NUP37", "PAK2", "RHOBTB2", "RHOV", "SPDL1")

# --- Load data ---
geneexp <- read.table("data/mRNA1.txt", header = TRUE, sep = "\t", 
                       check.names = FALSE, row.names = 1)
clinical <- read.table("data/clinical.txt", header = TRUE, sep = "\t")
clinical <- clinical[complete.cases(clinical[, 2]), ]
clinical <- clinical[, -c(2, 3)]

# Prepare expression data (tumor only, 12-char barcode)
colnames(geneexp) <- substr(colnames(geneexp), 1, 12)
geneexp <- geneexp[, !duplicated(colnames(geneexp))]

# Merge clinical + expression for signature genes
common <- intersect(clinical$id, colnames(geneexp))
clinical <- clinical[match(common, clinical$id), ]
geneexp_sig <- t(geneexp[signature_genes, common])

# Combine
data1 <- data.frame(
  survival_time = clinical$survival_time,
  status = clinical$status,
  geneexp_sig
)

# --- Multivariate Cox regression ---
gene_terms <- paste(signature_genes, collapse = " + ")
fmla <- as.formula(paste("Surv(survival_time, status) ~", gene_terms))
mycox <- coxph(fmla, data = data1)
coxResult <- summary(mycox)

# Extract results
multi_result <- data.frame(
  gene = names(coef(mycox)),
  coef = coef(mycox),
  HR = coxResult$conf.int[, "exp(coef)"],
  HR_95L = coxResult$conf.int[, "lower .95"],
  HR_95U = coxResult$conf.int[, "upper .95"],
  P = coxResult$coefficients[, "Pr(>|z|)"]
)

write.table(multi_result, "data/multivariate_cox_results.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Calculate risk score ---
# RiskScore = sum(coefficient_i * expression_i)
risk_score <- predict(mycox, type = "risk", newdata = data1)
risk_level <- ifelse(risk_score > median(risk_score), "High", "Low")

# Save risk scores
risk_df <- data.frame(
  id = common,
  risk_score = risk_score,
  risk_level = risk_level
)
write.table(risk_df, "data/risk_scores_final.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Median risk score cutoff:", median(risk_score), "\n")
cat("High risk:", sum(risk_level == "High"), 
    " Low risk:", sum(risk_level == "Low"), "\n")

# --- Forest plot of multivariate Cox ---
HRdf <- multi_result[, c("HR", "HR_95L", "HR_95U")]
rownames(HRdf) <- multi_result$gene

HR <- sprintf("%.4f", HRdf$HR)
HR_95L <- sprintf("%.4f", HRdf$HR_95L)
HR_95U <- sprintf("%.4f", HRdf$HR_95U)
P <- multi_result$P
P_fmt <- ifelse(P < 0.0001, "<0.0001", sprintf("%.4f", P))

wenben <- list(
  c(NA, multi_result$gene),
  append("P-value", P_fmt),
  append("Hazard ratio", paste0(HR, "(", HR_95L, "-", HR_95U, ")"))
)

pdf("figures/forest_multivariate_cox.pdf", width = 10, height = 8)
forestplot(wenben,
           rbind(rep(NA, 3), as.matrix(HRdf)),
           col = fpColors(box = "#ed4437", line = "black", summary = "black"),
           graphwidth = unit(50, "mm"),
           xlog = TRUE, lwd.ci = 2, boxsize = 0.3,
           xlab = "Hazard ratio")
dev.off()

# --- Identify independent prognostic factors ---
independent_genes <- multi_result[multi_result$P < 0.05, ]
cat("Independent prognostic factors (P < 0.05):\n")
print(independent_genes$gene)
