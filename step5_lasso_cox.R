###############################################################################
# Step 5: LASSO-Cox Regression (CORE MODEL BUILDING)

library(glmnet)
library(survival)
library(timeROC)  
library(ggplot2)

setwd("E:\\chendong RHO胰腺癌\\BJP修稿\\PAAD\\PAAD\\6")

# 1. 读取数据并预处理
mydata <- read.table("timeexpcd.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
mydata <- mydata[mydata$survival_time > 0, ]

v1 <- as.matrix(mydata[, c(3:ncol(mydata))])
v2 <- Surv(mydata$survival_time, mydata$status)

# ==============================================================================
# 1. 拟合 LASSO-Cox (统一采用 5 折交叉验证)
# ==============================================================================
set.seed(1234)

# 拟合 LASSO-Cox
myfit <- glmnet(v1, v2, family = "cox", alpha = 1)

# 5 折交叉验证 (与论文一致: nfolds = 5)
myfit2 <- cv.glmnet(v1, v2, family = "cox", nfolds = 5, alpha = 1)

# 导出 LASSO 路径与 5 折 CV 误差图
pdf("Figure3_LASSO_Path_and_CV.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(myfit, xvar = "lambda", label = TRUE)
abline(v = log(c(myfit2$lambda.min, myfit2$lambda.1se)), lty = "dashed", col = c("red", "blue"))
title("A. LASSO Coefficient Path", line = 2.5)

plot(myfit2)
abline(v = log(c(myfit2$lambda.min, myfit2$lambda.1se)), lty = "dashed", col = c("red", "blue"))
title("B. 5-Fold Cross-Validation Error", line = 2.5)
dev.off()
par(mfrow = c(1, 1))

# 计算 1 年期 (365天) 表观生存 AUC (Apparent 1-Year AUC)
pred_risk <- as.numeric(predict(myfit2, newx = v1, s = "lambda.min", type = "link"))
roc_orig  <- timeROC(T = mydata$survival_time, delta = mydata$status, 
                     marker = pred_risk, cause = 1, weighting = "marginal", 
                     times = 365, iid = FALSE)
obs_auc   <- roc_orig$AUC[2]

cat("\n================ 基础 LASSO-Cox 结果 ================\n")
cat("最佳 lambda.min :", myfit2$lambda.min, "\n")
cat("原始 1 年期表观 AUC (Apparent 1-Year AUC) :", round(obs_auc, 3), "\n")

# ==============================================================================
# 2. Bootstrap 重采样计算 1 年期 Optimism-Corrected AUC (B = 200)
# ==============================================================================
cat("\n[1/2] 正在运行 Bootstrap 评估 1 年期 AUC 乐观度 (B = 200)...\n")
B <- 200
optimism <- numeric(B)

set.seed(1234)
for (b in 1:B) {
  boot_idx <- sample(1:nrow(mydata), replace = TRUE)
  boot_data <- mydata[boot_idx, ]
  
  v1_boot <- v1[boot_idx, ]
  v2_boot <- v2[boot_idx]
  
  # 在 Bootstrap 集上独立进行 5 折 CV LASSO
  cv_boot <- try(cv.glmnet(v1_boot, v2_boot, family = "cox", alpha = 1, nfolds = 5), silent = TRUE)
  if (inherits(cv_boot, "try-error")) next
  
  # A. Bootstrap 模型在自举集上的 1 年期 AUC
  pred_boot_on_boot <- predict(cv_boot, newx = v1_boot, s = "lambda.min", type = "link")
  auc_boot <- timeROC(T = boot_data$survival_time, delta = boot_data$status, 
                      marker = pred_boot_on_boot, cause = 1, times = 365)$AUC[2]
  
  # B. Bootstrap 模型在原始数据集上的 1 年期 AUC
  pred_boot_on_orig <- predict(cv_boot, newx = v1, s = "lambda.min", type = "link")
  auc_orig <- timeROC(T = mydata$survival_time, delta = mydata$status, 
                      marker = pred_boot_on_orig, cause = 1, times = 365)$AUC[2]
  
  # 乐观度 = 自举 AUC - 原始 AUC
  optimism[b] <- auc_boot - auc_orig
}

mean_optimism <- mean(optimism, na.rm = TRUE)
corrected_auc <- obs_auc - mean_optimism

cat("Bootstrap 计算完成！\n")
cat("  - 1 年期平均乐观度 (Optimism):", round(mean_optimism, 3), "\n")
cat("  - 1 年期校正后 AUC (Optimism-Corrected AUC):", round(corrected_auc, 3), "\n")

# ==============================================================================
# 3. Permutation Test 置换检验 (N = 500, 经验 P 值加 1 修正)
# ==============================================================================
cat("\n[2/2] 正在运行 Permutation Test 验证信号非随机 (N = 500)...\n")
N_perm <- 500
perm_aucs <- numeric(N_perm)

set.seed(1234)
for (i in 1:N_perm) {
  perm_idx <- sample(1:nrow(mydata))
  perm_v2  <- v2[perm_idx]
  
  cv_perm <- try(cv.glmnet(v1, perm_v2, family = "cox", alpha = 1, nfolds = 5), silent = TRUE)
  if (inherits(cv_perm, "try-error")) {
    perm_aucs[i] <- 0.5
    next
  }
  
  pred_perm <- predict(cv_perm, newx = v1, s = "lambda.min", type = "link")
  roc_perm  <- timeROC(T = mydata$survival_time, delta = mydata$status[perm_idx], 
                       marker = pred_perm, cause = 1, times = 365)
  perm_aucs[i] <- roc_perm$AUC[2]
}

# 严格使用经验 P 值加 1 修正公式: P = (b + 1) / (N_perm + 1)
b_count <- sum(perm_aucs >= obs_auc, na.rm = TRUE)
p_val_corrected <- (b_count + 1) / (N_perm + 1)

# 导出 Permutation 图
perm_df <- data.frame(AUC = perm_aucs)
p_perm <- ggplot(perm_df, aes(x = AUC)) +
  geom_histogram(bins = 30, fill = "#7570B3", color = "white", alpha = 0.8) +
  geom_vline(xintercept = obs_auc, color = "firebrick", linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = obs_auc - 0.02, y = N_perm * 0.1, 
           label = paste0("Observed 1-Yr AUC = ", sprintf("%.3f", obs_auc), 
                          "\n(P = ", sprintf("%.3f", p_val_corrected), ")"), 
           color = "firebrick", hjust = 1, fontface = "bold") +
  labs(title = "Permutation Test Null Distribution (1-Year Survival)",
       subtitle = paste0("500 Permutations | Optimism-Corrected 1-Yr AUC = ", sprintf("%.3f", corrected_auc)),
       x = "Permuted 1-Year AUC", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        panel.grid.minor = element_blank())

ggsave("Figure3_Permutation_Test.pdf", plot = p_perm, width = 6.5, height = 4.5)
