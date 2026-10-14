###############################################################################
# Step 14: CIBERSORT Stacked Bar Plots
# Visualize immune cell composition across all samples
###############################################################################

library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
library(ggthemes)

# --- Load CIBERSORT results ---
cibersort_raw <- read.csv("data/CIBERSORT-Results.csv", row.names = 1, header = TRUE)

# --- Long-format data ---
dd1 <- cibersort_raw %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = 2:23,
               names_to = "CellType",
               values_to = "Composition")

plot.info <- dd1[, c("sample", "CellType", "Composition")]

# --- Full stacked bar plot ---
pdf("figures/immune_stacked_full.pdf", width = 15, height = 10)
ggbarplot(plot.info,
          x = "sample", y = "Composition",
          size = 0, fill = "CellType", color = "CellType") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 1),
    legend.position = "bottom"
  )
dev.off()

# --- Box plot of all cell types ---
pdf("figures/immune_boxplot_all.pdf", width = 15, height = 10)
ggboxplot(plot.info,
          x = "CellType", y = "Composition",
          color = "black", fill = "CellType",
          xlab = "", ylab = "Cell composition",
          main = "TME Cell Composition") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()
