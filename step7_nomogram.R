###############################################################################
# Step 7: Nomogram Construction + Calibration Curves
# Predict 1-, 2-, 3-year survival using risk score + tumor grade
###############################################################################

library(survival)
library(rms)

# --- Load data ---
riskcli <- read.table("data/riskcli_subset.txt", header = TRUE, sep = "\t")
# Expected columns: survival_time, status, age, sex, tumor_grade, 
#                   pathologic_stage, pathologic_T, pathologic_N, risk_score

# --- Data preparation ---
# Numeric conversion for categorical variables
conversion_rules <- list(
  sex = c("MALE" = 1, "FEMALE" = 2),
  tumor_grade = c("G1" = 1, "G2" = 2, "G3" = 3, "G4" = 4),
  pathologic_stage = c("Stage IA" = 1, "Stage IB" = 1, 
                       "Stage IIA" = 2, "Stage IIB" = 2, 
                       "Stage III" = 3, "Stage IV" = 4),
  pathologic_T = c("T1" = 1, "T2" = 2, "T3" = 3, "T4" = 4),
  pathologic_N = c("N0" = 1, "N1" = 2, "N1b" = 2)
)

# Apply conversions (simplified)
riskcli$grade_num <- as.numeric(factor(riskcli$tumor_grade, 
                                        levels = c("G1", "G2", "G3", "G4")))

# --- Survival function for nomogram ---
sur_1_year <- function(x) 0.1 * exp(-exp(x) * 1)   # placeholder
sur_2_year <- function(x) 0.1 * exp(-exp(x) * 2)
sur_3_year <- function(x) 0.1 * exp(-exp(x) * 3)

# Better approach: use rms baseline survival
dd <- datadist(riskcli)
options(datadist = "dd")

# Cox model for nomogram
cox <- cph(Surv(survival_time, status) ~ tumor_grade + risk_score,
           surv = TRUE, x = TRUE, y = TRUE, data = riskcli)

# --- Nomogram ---
pdf("figures/nomogram.pdf", width = 10, height = 7)
nom_sur <- nomogram(cox, 
                    fun = list(sur_1_year, sur_2_year, sur_3_year),
                    funlabel = c("1-year survival", "2-year survival", 
                                "3-year survival"),
                    maxscale = 100,
                    fun.at = seq(0.01, 0.99, by = 0.01))
plot(nom_sur)
dev.off()

# --- Calibration curves (Bootstrap B=1000) ---
# 1-year calibration
f1 <- cph(Surv(survival_time, status) ~ tumor_grade + risk_score,
          surv = TRUE, x = TRUE, y = TRUE, data = riskcli, time.inc = 365)

cal1 <- calibrate(f1, cmethod = "KM", method = "boot", 
                  u = 1, m = (nrow(riskcli) / 3), B = 1000)

# 2-year calibration
f2 <- cph(Surv(survival_time, status) ~ tumor_grade + risk_score,
          surv = TRUE, x = TRUE, y = TRUE, data = riskcli, time.inc = 730)

cal2 <- calibrate(f2, cmethod = "KM", method = "boot",
                  u = 2, m = (nrow(riskcli) / 3), B = 1000)

# 3-year calibration
f3 <- cph(Surv(survival_time, status) ~ tumor_grade + risk_score,
          surv = TRUE, x = TRUE, y = TRUE, data = riskcli, time.inc = 1095)

cal3 <- calibrate(f3, cmethod = "KM", method = "boot",
                  u = 3, m = (nrow(riskcli) / 3), B = 1000)

# --- Plot calibration curves ---
pdf("figures/calibration_curves.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
plot(cal1, xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "Nomogram-predicted 1-year survival",
     ylab = "Actual 1-year survival",
     main = "1-year calibration")
abline(0, 1, col = "gray", lty = 2)

plot(cal2, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Nomogram-predicted 2-year survival",
     ylab = "Actual 2-year survival",
     main = "2-year calibration")
abline(0, 1, col = "gray", lty = 2)

plot(cal3, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Nomogram-predicted 3-year survival",
     ylab = "Actual 3-year survival",
     main = "3-year calibration")
abline(0, 1, col = "gray", lty = 2)
dev.off()
