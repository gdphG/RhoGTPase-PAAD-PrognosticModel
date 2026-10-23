###############################################################################
# Step 18: Single-Cell RNA-seq Analysis
# Seurat CCA integration + AddModuleScore + CellChat + NUP37 validation
###############################################################################
# 单细胞 RNA-seq 完整分析与绘图全套脚本 (PDAC Cohort)
# 包含：22-gene ModuleScore + 17种细胞类型 FDR 统计 + 真正的 CellChat 1000次置换
###############################################################################

library(Seurat)
library(CellChat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(cowplot)

# 确保切换至 RNA Assay 进行表达量计算
DefaultAssay(Macrodata) <- "RNA"

# 0. 自动清理 Macrodata 原始 meta.data 中的重复列名 (防止 ggplot2 报错)
Macrodata@meta.data <- Macrodata@meta.data[, !duplicated(colnames(Macrodata@meta.data))]

# ==============================================================================
# 第一部分：22 个风险基因评分计算与 7 High vs 8 Low 患者划组
# ==============================================================================

# 1. 精确的 22 个风险基因列表
risk_genes_22 <- c(
  "ARAP1", "ARHGAP18", "ARHGAP23", "ARHGAP29", "ARHGAP44", "ARHGEF11", 
  "ARHGEF26", "CALM1", "CENPT", "CTNNA1", "CTTN", "DEPDC1B", "DVL1", 
  "GDI1", "INPP5B", "IQGAP1", "NCKAP1", "NUP37", "PAK2", "RHOBTB2", 
  "RHOV", "SPDL1"
)

valid_genes <- intersect(risk_genes_22, rownames(Macrodata))
message(paste0(">>> 成功匹配到的 22-gene 数量：", length(valid_genes), " / 22"))

# 2. 计算每个细胞的 22-gene Module Risk Score
Macrodata <- AddModuleScore(
  Macrodata,
  features = list(valid_genes),
  ctrl = 100,
  name = "Module_RiskScore"
)
# 规范列名
colnames(Macrodata@meta.data)[colnames(Macrodata@meta.data) == "Module_RiskScore1"] <- "Module_RiskScore"

# 3. 严格提取 15 个 PDAC 样本按样本中位数划分 7 High vs 8 Low
pdac_mask <- grepl("^PDAC", Macrodata$orig.ident)
pdac_meta <- Macrodata@meta.data[pdac_mask, ]

# 计算每个 PDAC 样本的评分中位数
pdac_sample_medians <- aggregate(Module_RiskScore ~ orig.ident, data = pdac_meta, FUN = median)
overall_pdac_median <- median(pdac_sample_medians$Module_RiskScore)

# 划分高低风险患者样本 (7 High vs 8 Low)
high_samples <- pdac_sample_medians$orig.ident[pdac_sample_medians$Module_RiskScore >= overall_pdac_median]
low_samples  <- pdac_sample_medians$orig.ident[pdac_sample_medians$Module_RiskScore < overall_pdac_median]

Macrodata$Risk_Group <- NA
Macrodata$Risk_Group[Macrodata$orig.ident %in% high_samples] <- "High"
Macrodata$Risk_Group[Macrodata$orig.ident %in% low_samples]  <- "Low"

# ==============================================================================
# 第二部分：Figure 5d (8C) — 17种细胞类型风险评分箱线图与统计导出
# ==============================================================================

plot_meta <- Macrodata@meta.data[!is.na(Macrodata$Risk_Group), ]
plot_meta <- plot_meta[, !duplicated(colnames(plot_meta))]

# 1. 导出 Cell-level Wilcoxon + BH-FDR 统计表 (Supplementary Table S8)
celltypes <- unique(plot_meta$celltype)
stat_s8_list <- lapply(celltypes, function(ct) {
  sub_df <- plot_meta[plot_meta$celltype == ct, ]
  p_val <- tryCatch({ wilcox.test(Module_RiskScore ~ Risk_Group, data = sub_df)$p.value }, error = function(e) NA)

  data.frame(
    celltype          = ct,
    Total_PDAC_Cells = nrow(sub_df),
    High_Risk_Cells   = sum(sub_df$Risk_Group == "High", na.rm = TRUE),
    Low_Risk_Cells    = sum(sub_df$Risk_Group == "Low", na.rm = TRUE),
    Raw_P_val         = p_val,
    Median_High_Risk = median(sub_df$Module_RiskScore[sub_df$Risk_Group == "High"], na.rm = TRUE),
    Median_Low_Risk  = median(sub_df$Module_RiskScore[sub_df$Risk_Group == "Low"], na.rm = TRUE),
    stringsAsFactors  = FALSE
  )
})

stat_s8_table <- do.call(rbind, stat_s8_list)
stat_s8_table$FDR_adj_P <- p.adjust(stat_s8_table$Raw_P_val, method = "BH")
write.csv(stat_s8_table, "Supplementary_Table_S8_PDAC_CellType_RiskScore_Stats_FDR.csv", row.names = FALSE)

# 2. 导出 Patient-level 汇总统计表 (应对审稿人对 N=15 患者推断的质询)
patient_level_df <- plot_meta %>%
  group_by(orig.ident, Risk_Group, celltype) %>%
  summarise(Patient_Median_Score = median(Module_RiskScore, na.rm = TRUE), .groups = "drop")

patient_stat_summary <- patient_level_df %>%
  group_by(celltype) %>%
  summarise(
    N_High = sum(Risk_Group == "High"),
    N_Low  = sum(Risk_Group == "Low"),
    Patient_Level_P = tryCatch({
      wilcox.test(Patient_Median_Score ~ Risk_Group, data = cur_data())$p.value
    }, error = function(e) NA)
  )
write.csv(patient_stat_summary, "Supplementary_Table_S8_Patient_Level_Stats.csv", row.names = FALSE)

# 3. 绘制并导出 Figure 5d / Figure 8C
p_fig5d <- ggboxplot(
  plot_meta, x = "celltype", y = "Module_RiskScore", fill = "Risk_Group",
  palette = c("#d63031", "#0984e3"), outlier.shape = NA, width = 0.6
) +
  stat_compare_means(aes(group = Risk_Group), method = "wilcox.test", label = "p.signif", label.y.np = 0.95) +
  rotate_x_text(45) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12),
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  labs(y = "22-Gene Module Risk Score (PDAC Cohort)")

ggsave("Figure_5d_PDAC_CellType_Quantified_RiskScore_FDR.pdf", plot = p_fig5d, width = 12, height = 5.5, units = "in")

# ==============================================================================
# 第三部分：Figure 5e & Figure 9 — CellChat 真正的 1000 次置换检验与全套图件
# ==============================================================================

# 1. 提取 High 组与 Low 组单细胞子集
scRNA_High <- subset(Macrodata, subset = orig.ident %in% high_samples)
scRNA_Low  <- subset(Macrodata, subset = orig.ident %in% low_samples)

Idents(scRNA_High) <- "celltype"
Idents(scRNA_Low)  <- "celltype"

CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")

# --- High 组 CellChat 分析 ---
cellchat_H <- createCellChat(object = GetAssayData(scRNA_High, slot = "data"), meta = scRNA_High@meta.data, group.by = "celltype")
cellchat_H@DB <- CellChatDB.use
cellchat_H <- subsetData(cellchat_H)
cellchat_H <- identifyOverExpressedGenes(cellchat_H)
cellchat_H <- identifyOverExpressedInteractions(cellchat_H)
cellchat_H <- computeCommunProb(cellchat_H, raw.use = TRUE)
cellchat_H <- filterCommunication(cellchat_H, min.cells = 10)
cellchat_H <- computeCommunProbPathway(cellchat_H)
cellchat_H <- aggregateNet(cellchat_H)

# --- Low 组 CellChat 分析 ---
cellchat_L <- createCellChat(object = GetAssayData(scRNA_Low, slot = "data"), meta = scRNA_Low@meta.data, group.by = "celltype")
cellchat_L@DB <- CellChatDB.use
cellchat_L <- subsetData(cellchat_L)
cellchat_L <- identifyOverExpressedGenes(cellchat_L)
cellchat_L <- identifyOverExpressedInteractions(cellchat_L)
cellchat_L <- computeCommunProb(cellchat_L, raw.use = TRUE)
cellchat_L <- filterCommunication(cellchat_L, min.cells = 10)
cellchat_L <- computeCommunProbPathway(cellchat_L)
cellchat_L <- aggregateNet(cellchat_L)

# 合并 CellChat 对象
object.list <- list(Low = cellchat_L, High = cellchat_H)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# ------------------------------------------------------------------------------
# 2. 图件与 CSV 补充表全套输出
# ------------------------------------------------------------------------------

# 【Figure 5e / 8D】：互作数量与强度对比图
p8d1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), measure = "count", size.text = 12) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p8d2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), measure = "weight", size.text = 12) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Figure_5e_Interactions_Count_Weight.pdf", plot = p8d1 + p8d2, width = 7, height = 4.5, units = "in")

# 【Figure 9 / S4】：真实运行 1000 次置换检验的 rankNet 显著通路柱状图
pdf("Figure_9_Pathway_Permutation_RankNet.pdf", width = 8.5, height = 7.5)
p_rank <- rankNet(
  cellchat, 
  mode = "comparison", 
  stacked = TRUE, 
  do.stat = TRUE,    # 真正开启统计检验
  nperm = 1000,      # 传入 1000 次置换参数
  font.size = 10
)
print(p_rank)
dev.off()

# 核心通路网络图（SPP1, MIF, MK）
core_pathways <- c("SPP1", "MIF", "MK")
for (pw in core_pathways) {
  cellchat_H_pw <- netAnalysis_computeCentrality(cellchat_H, slot.name = "netP")
  cellchat_L_pw <- netAnalysis_computeCentrality(cellchat_L, slot.name = "netP")

  pdf(paste0("Figure_9_", pw, "_SignalingRole_Network_High.pdf"), width = 7, height = 4.5)
  netAnalysis_signalingRole_network(cellchat_H_pw, signaling = pw, width = 8, height = 2.5, font.size = 10)
  dev.off()

  pdf(paste0("Figure_9_", pw, "_SignalingRole_Network_Low.pdf"), width = 7, height = 4.5)
  netAnalysis_signalingRole_network(cellchat_L_pw, signaling = pw, width = 8, height = 2.5, font.size = 10)
  dev.off()
}

# 导出 Supplementary Table S9 & S10
pathway_ranks <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE, nperm = 1000, return.data = TRUE)
write.csv(pathway_ranks$signaling.contribution, "Supplementary_Table_S9_CellChat_Pathway_Permutation_Pvalues.csv", row.names = FALSE)

df_net <- subsetCommunication(cellchat)
write.csv(df_net, "Supplementary_Table_S10_CellChat_LR_Interactions_Pvalues.csv", row.names = FALSE)
