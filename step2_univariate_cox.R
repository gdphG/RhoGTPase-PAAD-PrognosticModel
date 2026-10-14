###############################################################################
# Step 2: Univariate Cox Regression Analysis
# Screen prognosis-associated genes from Rho GTPase gene list
# BH correction for P-values, threshold: adjusted P < 0.05
###############################################################################

library(survival)
library(forestplot)

# --- Load expression data and clinical data ---
geneexp <- read.table("data/mRNA1.txt", header = TRUE, sep = "\t", 
                       check.names = FALSE, row.names = 1)

clinical <- read.table("data/clinical.txt", header = TRUE, sep = "\t")
# Keep only samples with complete survival data
clinical <- clinical[complete.cases(clinical[, 2]), ]
clinical <- clinical[, -c(2, 3)]  # Remove non-essential columns

# --- Extract tumor-only samples (TCGA barcode) ---
colnames(geneexp) <- substr(colnames(geneexp), 1, 12)
geneexp <- geneexp[, !duplicated(colnames(geneexp))]

# Merge with clinical data
common_samples <- intersect(colnames(geneexp), clinical$id)
geneexp <- geneexp[, common_samples]
clinical <- clinical[match(common_samples, clinical$id), ]

# --- Univariate Cox regression ---
# Define Rho GTPase gene list (124 genes from literature)
rho_genes <- read.table("data/rho_gtpase_genes.txt", header = FALSE)[, 1]
# Alternative: use all genes in expression matrix for screening

# Cox regression function
cox_univariate <- function(gene_name) {
  expr <- as.numeric(geneexp[gene_name, ])
  df <- data.frame(
    survival_time = clinical$survival_time,
    status = clinical$status,
    expr = expr
  )
  fmla <- as.formula(Surv(survival_time, status) ~ expr)
  mycox <- coxph(fmla, data = df)
  coxResult <- summary(mycox)
  
  data.frame(
    gene = gene_name,
    HR = coxResult$conf.int[, "exp(coef)"],
    HR_95L = coxResult$conf.int[, "lower .95"],
    HR_95U = coxResult$conf.int[, "upper .95"],
    P = coxResult$coefficients[, "Pr(>|z|)"]
  )
}

# Run for all genes
results <- do.call(rbind, lapply(rho_genes, cox_univariate))

# --- BH correction ---
results$adjusted_P <- p.adjust(results$P, method = "BH")

# Filter significant genes
sig_genes <- results[results$adjusted_P < 0.05, ]
write.table(sig_genes, "data/univariate_cox_significant.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Forest plot ---
HRdf <- as.matrix(sig_genes[, c("HR", "HR_95L", "HR_95U")])
rownames(HRdf) <- sig_genes$gene

HR <- sprintf("%.4f", as.numeric(HRdf[, "HR"]))
HR_95L <- sprintf("%.4f", as.numeric(HRdf[, "HR_95L"]))
HR_95U <- sprintf("%.4f", as.numeric(HRdf[, "HR_95U"]))
P <- sig_genes$adjusted_P
P_fmt <- ifelse(P < 0.0001, "<0.0001", sprintf("%.4f", P))

wenben <- list(
  c(NA, rownames(HRdf)),
  append("P", P_fmt),
  append("Hazard ratio", paste0(HR, "(", HR_95L, "-", HR_95U, ")"))
)

pdf("figures/forest_univariate.pdf", width = 10, height = 8)
forestplot(wenben,
           rbind(rep(NA, 3), HRdf),
           col = fpColors(box = "#4974a4", line = "black", summary = "black"),
           graphwidth = unit(50, "mm"),
           xlog = TRUE, lwd.ci = 2, boxsize = 0.3,
           xlab = "Hazard ratio")
dev.off()

cat("Number of significant prognosis-associated genes:", nrow(sig_genes), "\n")
