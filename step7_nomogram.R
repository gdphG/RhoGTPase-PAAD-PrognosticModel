###############################################################################
# Step 7: Nomogram Construction + Calibration Curves
# Predict 1-, 2-, 3-year survival using risk score + tumor grade
###############################################################################

library(survival)
library(rms)

# --- 1. Load and prepare data ---
riskcli <- read.table("data/riskcli_subset.txt", header = TRUE, sep = "\t")

# Ensure survival_time is numeric and properly formatted
riskcli$survival_time <- as.numeric(riskcli$survival_time)
riskcli$status <- as.numeric(riskcli$status)

# Convert factor variables
riskcli$tumor_grade <- factor(riskcli$tumor_grade, levels = c("G1", "G2", "G3", "G4"))

# --- 2. rms Data Environment Setup ---
dd <- datadist(riskcli)
options(datadist = "dd")

# --- 3. Fit cph model for Nomogram ---
# Time unit is in DAYS: 1-Year = 365, 2-Year = 730, 3-Year = 1095
cox <- cph(Surv(survival_time, status) ~ tumor_grade + risk_score,
           surv = TRUE, x = TRUE, y = TRUE, data = riskcli)

# --- 4. Dynamically Bind Standard Survival Functions (Fixes 0.1 Cap Error) ---
# Use the rms Survival() factory function to evaluate S(t) directly from the model
surv_obj <- Survival(cox)
sur_1_year <- function(x) surv_obj(365, x)   # 1-Year (365 Days)
sur_2_year <- function(x) surv_obj(730, x)   # 2-Year (730 Days)
sur_3_year <- function(x) surv_obj(1095, x)  # 3-Year (1095 Days)

# --- 5. Generate Nomogram ---
pdf("figures/nomogram.pdf", width = 10, height = 7)
nom_sur <- nomogram(cox, 
                    fun = list(sur_1_year, sur_2_year, sur_3_year),
                    funlabel = c("1-Year Survival Probability", 
                                 "2-Year Survival Probability", 
                                 "3-Year Survival Probability"),
                    maxscale = 100,
                    fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9)) # Clear tick marks from 0.1 to 0.9
plot(nom_sur)
dev.off()

# --- 6. Calibration Curves (Harmonized Time Unit: u = days) ---
# 1-Year Calibration (u = 365 days)
f1 <- cph(Surv(survival_time, status) ~ tumor_grade + risk_score,
          surv = TRUE, x = TRUE, y = TRUE, data = riskcli, time.inc = 365)
cal1 <- calibrate(f1, cmethod = "KM", method = "boot", 
                  u = 365, m = (nrow(riskcli) / 3), B = 1000)

# 2-Year Calibration (u = 730 days)
f2 <- cph(Surv(survival_time, status) ~ tumor_grade + risk_score,
          surv = TRUE, x = TRUE, y = TRUE, data = riskcli, time.inc = 730)
cal2 <- calibrate(f2, cmethod = "KM", method = "boot", 
                  u = 730, m = (nrow(riskcli) / 3), B = 1000)

# 3-Year Calibration (u = 1095 days)
f3 <- cph(Surv(survival_time, status) ~ tumor_grade + risk_score,
          surv = TRUE, x = TRUE, y = TRUE, data = riskcli, time.inc = 1095)
cal3 <- calibrate(f3, cmethod = "KM", method = "boot", 
                  u = 1095, m = (nrow(riskcli) / 3), B = 1000)

# --- 7. Plot Calibration Curves ---
pdf("figures/calibration_curves.pdf", width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

plot(cal1, xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "Nomogram-Predicted 1-Year Survival",
     ylab = "Actual 1-Year Survival",
     main = "1-Year Calibration")

plot(cal2, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Nomogram-Predicted 2-Year Survival",
     ylab = "Actual 2-Year Survival",
     main = "2-Year Calibration")

plot(cal3, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Nomogram-Predicted 3-Year Survival",
     ylab = "Actual 3-Year Survival",
     main = "3-Year Calibration")

dev.off()
