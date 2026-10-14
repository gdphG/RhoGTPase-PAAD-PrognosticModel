###############################################################################
# Step 9: Subgroup Cox Analysis + Clinical Feature Forest Plots
# Univariate and multivariate Cox with clinical variables
###############################################################################

library(survival)
library(forestplot)
library(dplyr)

# --- Load and prepare data ---
riskcli <- read.table("data/riskcli_subset.txt", header = TRUE, sep = "\t")

# Numeric conversion for clinical variables
conversion_rules <- list(
  sex = c("MALE" = 1, "FEMALE" = 2),
  tumor_grade = c("G1" = 1, "G2" = 2, "G3" = 3, "G4" = 4),
  pathologic_stage = c("Stage IA" = 1, "Stage IB" = 1, 
                       "Stage IIA" = 2, "Stage IIB" = 2, 
                       "Stage III" = 3, "Stage IV" = 4),
  pathologic_T = c("T1" = 1, "T2" = 2, "T3" = 3, "T4" = 4),
  pathologic_N = c("N0" = 1, "N1" = 2, "N1b" = 2),
  risk_level = c("High" = 2, "Low" = 1)
)

riskcli <- riskcli %>%
  mutate(
    sex_numeric = case_when(
      sex %in% names(conversion_rules$sex) ~ as.integer(conversion_rules$sex[sex]),
      TRUE ~ NA_integer_
    ),
    tumor_grade_numeric = case_when(
      tumor_grade %in% names(conversion_rules$tumor_grade) ~ as.integer(conversion_rules$tumor_grade[tumor_grade]),
      TRUE ~ NA_integer_
    ),
    pathologic_stage_numeric = case_when(
      pathologic_stage %in% names(conversion_rules$pathologic_stage) ~ as.integer(conversion_rules$pathologic_stage[pathologic_stage]),
      TRUE ~ NA_integer_
    ),
    pathologic_T_numeric = case_when(
      pathologic_T %in% names(conversion_rules$pathologic_T) ~ as.integer(conversion_rules$pathologic_T[pathologic_T]),
      TRUE ~ NA_integer_
    ),
    pathologic_N_numeric = case_when(
      pathologic_N %in% names(conversion_rules$pathologic_N) ~ as.integer(conversion_rules$pathologic_N[pathologic_N]),
      TRUE ~ NA_integer_
    ),
    risk_level_numeric = case_when(
      risk_level %in% names(conversion_rules$risk_level) ~ as.integer(conversion_rules$risk_level[risk_level]),
      TRUE ~ NA_integer_
    )
  )

row.names(riskcli) <- riskcli[, 1]

# --- Univariate Cox for each clinical variable ---
data1 <- riskcli
coxf <- function(x) {
  fmla1 <- as.formula(paste0("Surv(survival_time, status) ~ data1[, '", x, "']"))
  mycox <- coxph(fmla1, data = data1)
}

newdf <- data.frame()
for (a in colnames(data1[, 3:ncol(data1)])) {
  mycox <- coxf(a)
  coxResult <- summary(mycox)
  newdf <- rbind(newdf, cbind(
    id = a,
    HR = coxResult$conf.int[, "exp(coef)"],
    HR_95L = coxResult$conf.int[, "lower .95"],
    HR_95U = coxResult$conf.int[, "upper .95"],
    P = coxResult$coefficients[, "Pr(>|z|)"]
  ))
}
write.table(newdf, "data/subgroup_univariate_cox.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)

# --- Forest plot (univariate) ---
HRdf <- as.matrix(newdf[, -1])
rownames(HRdf) <- newdf$id
HR <- sprintf("%.4f", as.numeric(HRdf[, "HR"]))
P <- as.numeric(HRdf[, "P"])
P_fmt <- ifelse(P < 0.0001, "<0.0001", sprintf("%.4f", P))

wenben <- list(
  c(NA, rownames(HRdf)),
  append("P", P_fmt),
  append("Hazard ratio", paste0(HR, "(", HRdf[, "HR_95L"], "-", HRdf[, "HR_95U"], ")"))
)

pdf("figures/forest_clinical_univariate.pdf", width = 8, height = 6)
forestplot(wenben,
           rbind(rep(NA, 3), HRdf[, 1:3]),
           col = fpColors(box = "#4974a4", line = "black", summary = "black"),
           graphwidth = unit(50, "mm"),
           xlog = TRUE, lwd.ci = 2, boxsize = 0.3,
           xlab = "Hazard ratio")
dev.off()

# --- Multivariate Cox ---
fmla2 <- as.formula(Surv(survival_time, status) ~ Age + Sex + Grade + Stage + T + N + Riskscore)
mycox2 <- coxph(fmla2, data = data1)
coxResult2 <- summary(mycox2)

multi_df <- cbind(
  id = row.names(coxResult2$conf.int),
  HR = coxResult2$conf.int[, "exp(coef)"],
  HR_95L = coxResult2$conf.int[, "lower .95"],
  HR_95U = coxResult2$conf.int[, "upper .95"],
  P = coxResult2$coefficients[, "Pr(>|z|)"]
)
write.table(multi_df, "data/subgroup_multivariate_cox.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

# --- Forest plot (multivariate) ---
HRdf2 <- as.matrix(multi_df[, -1])
rownames(HRdf2) <- multi_df$id
P2 <- as.numeric(HRdf2[, "P"])
P_fmt2 <- ifelse(P2 < 0.001, "<0.001", sprintf("%.4f", P2))

wenben2 <- list(
  c(NA, rownames(HRdf2)),
  append("P-value", P_fmt2),
  append("Hazard ratio", paste0(HRdf2[, "HR"], "(", HRdf2[, "HR_95L"], "-", HRdf2[, "HR_95U"], ")"))
)

pdf("figures/forest_clinical_multivariate.pdf", width = 10, height = 6)
forestplot(wenben2,
           rbind(rep(NA, 3), HRdf2[, 1:3]),
           col = fpColors(box = "#ed4437", line = "black", summary = "black"),
           graphwidth = unit(50, "mm"),
           xlog = TRUE, lwd.ci = 2, boxsize = 0.3,
           xlab = "Hazard ratio")
dev.off()
