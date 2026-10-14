###############################################################################
# Step 11: GSEA (Gene Set Enrichment Analysis)
# KEGG and GO enrichment between High vs. Low risk groups
###############################################################################

library(limma)
library(clusterProfiler)
library(enrichplot)
library(dplyr)

# --- Load data ---
exp_data <- read.table("data/mRNA1.txt", header = TRUE, sep = "\t", 
                        row.names = 1, check.names = FALSE)
risk_data <- read.table("data/risk_scores_final.txt", header = TRUE, sep = "\t")
row.names(risk_data) <- risk_data$id

# Merge and prepare
merged <- merge(t(exp_data), risk_data, by = "row.names")
row.names(merged) <- merged$Row.names
merged <- merged[, -1]

# --- Prepare expression matrix (tumor only) ---
group <- sapply(strsplit(colnames(exp_data), "-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
group <- gsub("2", "1", group)
exp_tumor <- exp_data[, group == 1]
exp_tumor <- t(avereps(t(exp_tumor)))

# Align samples
common_samples <- intersect(colnames(exp_tumor), rownames(merged))
exp_tumor <- exp_tumor[, common_samples]

# --- Calculate logFC (High vs. Low) ---
dataL <- exp_tumor[, rownames(merged[merged$risk_level == "Low", ])]
dataH <- exp_tumor[, rownames(merged[merged$risk_level == "High", ])]
meanL <- rowMeans(dataL)
meanH <- rowMeans(dataH)
meanL[meanL < 0.00001] <- 0.00001
meanH[meanH < 0.00001] <- 0.00001
logFC <- log2(meanH) - log2(meanL)
logFC <- sort(logFC, decreasing = TRUE)

# ==========================================================================
# KEGG GSEA
# ==========================================================================
gmt_kegg <- read.gmt("data/c2.cp.kegg.v7.4.symbols.gmt")
kk_kegg <- GSEA(logFC, TERM2GENE = gmt_kegg, pvalueCutoff = 1)
kkTab_kegg <- as.data.frame(kk_kegg)
kkTab_kegg <- kkTab_kegg[kkTab_kegg$pvalue < 0.05, ]
write.table(kkTab_kegg, "data/GSEA_KEGG_results.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE)

# Plot top 5 enriched pathways for High risk
pdf("figures/GSEA_KEGG_high_risk.pdf", width = 7, height = 5.5)
if (nrow(kkTab_kegg[kkTab_kegg$NES > 0, ]) >= 5) {
  showTerm <- rownames(kkTab_kegg[kkTab_kegg$NES > 0, ])[1:5]
  print(gseaplot2(kk_kegg, showTerm, base_size = 8, title = "Enriched in high risk group"))
}
dev.off()

# Plot for Low risk
pdf("figures/GSEA_KEGG_low_risk.pdf", width = 7, height = 5.5)
if (nrow(kkTab_kegg[kkTab_kegg$NES < 0, ]) >= 5) {
  showTerm <- rownames(kkTab_kegg[kkTab_kegg$NES < 0, ])[1:5]
  print(gseaplot2(kk_kegg, showTerm, base_size = 8, title = "Enriched in low risk group"))
}
dev.off()

# ==========================================================================
# GO GSEA
# ==========================================================================
gmt_go <- read.gmt("data/c5.go.v7.4.symbols.gmt")
kk_go <- GSEA(logFC, TERM2GENE = gmt_go, pvalueCutoff = 1)
kkTab_go <- as.data.frame(kk_go)
kkTab_go <- kkTab_go[kkTab_go$pvalue < 0.05, ]
write.table(kkTab_go, "data/GSEA_GO_results.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE)

# Plot
pdf("figures/GSEA_GO_high_risk.pdf", width = 7, height = 5.5)
if (nrow(kkTab_go[kkTab_go$NES > 0, ]) >= 5) {
  showTerm <- rownames(kkTab_go[kkTab_go$NES > 0, ])[1:5]
  print(gseaplot2(kk_go, showTerm, base_size = 8, title = "Enriched in high risk group"))
}
dev.off()

pdf("figures/GSEA_GO_low_risk.pdf", width = 7, height = 5.5)
if (nrow(kkTab_go[kkTab_go$NES < 0, ]) >= 5) {
  showTerm <- rownames(kkTab_go[kkTab_go$NES < 0, ])[1:5]
  print(gseaplot2(kk_go, showTerm, base_size = 8, title = "Enriched in low risk group"))
}
dev.off()
