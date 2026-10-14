###############################################################################
# Step 15: Somatic Mutation Landscape Analysis
# Using maftools + TCGAmutations for PAAD mutation data
###############################################################################

library(tidyverse)
library(maftools)
library(TCGAmutations)

# --- Download and load PAAD mutation data ---
TCGAmutations::tcga_available()
mut_data <- TCGAmutations::tcga_load(study = "PAAD")

# --- Prepare mutation matrix ---
a <- mut_data@data
a <- a[, c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")] %>%
  mutate(Tumor_Sample_Barcode = substring(Tumor_Sample_Barcode, 1, 16))

# Count mutations per gene
gene <- unique(a$Hugo_Symbol)
sample <- unique(a$Tumor_Sample_Barcode)

mat <- as.data.frame(matrix("", length(gene), length(sample),
                             dimnames = list(gene, sample)))
mat_0_1 <- as.data.frame(matrix(0, length(gene), length(sample),
                                 dimnames = list(gene, sample)))

for (i in 1:nrow(a)) {
  mat[a[i, 1], a[i, 3]] <- a[i, 2]
  mat_0_1[a[i, 1], a[i, 3]] <- 1
}

gene_count <- data.frame(gene = rownames(mat_0_1),
                          count = rowSums(mat_0_1)) %>%
  arrange(desc(count))

save(mat, mat_0_1, gene_count, file = "data/PAAD_mutation_data.rda")

# --- Oncoplot (overall top 10) ---
pdf("figures/mutation_oncoplot_top10.pdf", width = 10, height = 8)
oncoplot(maf = mut_data, top = 10, fontSize = 0.6, 
         showTumorSampleBarcodes = FALSE)
dev.off()

# --- Split by risk group ---
group <- read.table("data/risk_scores_final.txt", header = TRUE, sep = "\t")
group <- group[, c("id", "risk_level")]
colnames(group) <- c("Tumor_Sample_Barcode", "group")
a$Tumor_Sample_Barcode <- substring(a$Tumor_Sample_Barcode, 1, 12)

maf_high <- a[a$Tumor_Sample_Barcode %in% group$Tumor_Sample_Barcode[group$group == "High"], ]
maf_low <- a[a$Tumor_Sample_Barcode %in% group$Tumor_Sample_Barcode[group$group == "Low"], ]

maf_obj_high <- read.maf(maf = maf_high)
maf_obj_low <- read.maf(maf = maf_low)

# --- Oncoplots by risk group ---
pdf("figures/mutation_oncoplot_high_risk.pdf", width = 10, height = 8)
oncoplot(maf = maf_obj_high, top = 10, fontSize = 0.6,
         showTumorSampleBarcodes = FALSE)
dev.off()

pdf("figures/mutation_oncoplot_low_risk.pdf", width = 10, height = 8)
oncoplot(maf = maf_obj_low, top = 20, fontSize = 0.6,
         showTumorSampleBarcodes = FALSE)
dev.off()

# --- TMB comparison ---
tmb_high <- length(which(mat_0_1[, group$Tumor_Sample_Barcode[group$group == "High"]] > 0))
tmb_low <- length(which(mat_0_1[, group$Tumor_Sample_Barcode[group$group == "Low"]] > 0))

# Wilcoxon test
tmb_by_sample <- colSums(mat_0_1)
names(tmb_by_sample) <- colnames(mat_0_1)
tmb_df <- data.frame(
  sample = names(tmb_by_sample),
  TMB = tmb_by_sample
)
tmb_df$group <- group$group[match(tmb_df$sample, group$Tumor_Sample_Barcode)]
tmb_df <- tmb_df[!is.na(tmb_df$group), ]

wilcox.test(TMB ~ group, data = tmb_df)

pdf("figures/TMB_comparison.pdf", width = 5, height = 5)
ggplot(tmb_df, aes(x = group, y = TMB, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  theme_bw() +
  ggtitle("Tumor Mutation Burden: High vs. Low Risk")
dev.off()
