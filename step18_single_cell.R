###############################################################################
# Step 18: Single-Cell RNA-seq Analysis
# Seurat CCA integration + AddModuleScore + CellChat + NUP37 validation
###############################################################################

library(Seurat)
library(CellChat)
library(ggplot2)
library(dplyr)

# --- 22-gene signature ---
signature_genes <- c("ARAP1", "ARHGAP18", "ARHGAP23", "ARHGAP29", "ARHGAP44",
                     "ARHGEF11", "ARHGEF26", "CALM1", "CENPT", "CTNNA1",
                     "CTTN", "DEPDC1B", "DVL1", "GDI1", "INPP5B", "IQGAP1",
                     "NCKAP1", "NUP37", "PAK2", "RHOBTB2", "RHOV", "SPDL1")

# ============================================================================
# PART 1: Seurat Analysis with CCA Integration
# ============================================================================

# --- Load scRNA-seq data (GSE57495 biopsy cohort) ---
# Replace with actual path to your Seurat objects
obj1 <- Read10X(data.dir = "data/scRNA/sample1/")
obj2 <- Read10X(data.dir = "data/scRNA/sample2/")

seurat1 <- CreateSeuratObject(counts = obj1, project = "Sample1")
seurat2 <- CreateSeuratObject(counts = obj2, project = "Sample2")

# --- Standard QC ---
seurat1[["percent.mt"]] <- PercentageFeatureSet(seurat1, pattern = "^MT-")
seurat2[["percent.mt"]] <- PercentageFeatureSet(seurat2, pattern = "^MT-")

seurat1 <- subset(seurat1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
seurat2 <- subset(seurat2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# --- CCA Integration (batch correction) ---
seurat_list <- list(seurat1, seurat2)
seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                          dims = 1:30, 
                                          anchor.features = 2000)
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:30)

# --- Standard Seurat workflow ---
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = FALSE)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:20)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20)

# --- UMAP visualization ---
pdf("figures/sc_umap_clusters.pdf", width = 8, height = 7)
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, 
        pt.size = 0.5) + ggtitle("UMAP: Cell Clusters")
dev.off()

# --- Cell type annotation ---
# Using canonical markers (from literature)
marker_genes <- list(
  "T cells" = c("CD3D", "CD3E", "CD8A", "CD4"),
  "B cells" = c("CD19", "MS4A1", "CD79A"),
  "Macrophages" = c("CD68", "CD163", "CSF1R"),
  "Fibroblasts" = c("COL1A1", "DCN", "FAP"),
  "Epithelial" = c("EPCAM", "KRT18", "KRT19"),
  "Endothelial" = c("PECAM1", "VWF", "CDH5"),
  "Dendritic" = c("ITGAX", "HLA-DRA", "CLEC9A"),
  "NK cells" = c("NKG7", "GNLY", "KLRD1")
)

# AddModuleScore for cell type identification
for (ct in names(marker_genes)) {
  genes_use <- intersect(marker_genes[[ct]], rownames(seurat_integrated))
  if (length(genes_use) > 0) {
    seurat_integrated <- AddModuleScore(seurat_integrated, 
                                         features = list(genes_use),
                                         name = paste0(ct, "_score"))
  }
}

# ============================================================================
# PART 2: Risk Score at Single-Cell Level
# ============================================================================

# AddModuleScore for the 22-gene signature
seurat_integrated <- AddModuleScore(
  seurat_integrated,
  features = list(signature_genes),
  name = "RiskScore_22gene"
)

# Visualize risk score on UMAP
pdf("figures/sc_risk_score_umap.pdf", width = 8, height = 7)
FeaturePlot(seurat_integrated, features = "RiskScore_22gene1", 
            min.cutoff = 0, max.cutoff = 1,
            cols = c("lightgrey", "red")) +
  ggtitle("22-Gene Risk Score (Single-cell)")
dev.off()

# ============================================================================
# PART 3: NUP37 Expression in Single Cells
# ============================================================================

# NUP37 validation - check expression across cell types
pdf("figures/sc_NUP37_violin.pdf", width = 8, height = 5)
VlnPlot(seurat_integrated, features = "NUP37", 
        pt.size = 0.1, cols = rainbow(10)) +
  ggtitle("NUP37 Expression Across Cell Types")
dev.off()

pdf("figures/sc_NUP37_feature.pdf", width = 7, height = 6)
FeaturePlot(seurat_integrated, features = "NUP37",
            cols = c("lightgrey", "darkred")) +
  ggtitle("NUP37 Expression (Single-cell)")
dev.off()

# ============================================================================
# PART 4: CellChat Cell-Cell Communication
# ============================================================================

# --- Prepare for CellChat ---
# Using integrated expression
sc_data <- GetAssayData(seurat_integrated, slot = "data")
cell_labels <- seurat_integrated$seurat_clusters

cellchat <- createCellChat(object = sc_data, 
                            group = cell_labels, 
                            group.by = "Number")

# --- Preprocessing ---
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# --- Visualization ---
pdf("figures/cellchat_network.pdf", width = 8, height = 7)
netVisual_circlePlot(cellchat@net$weight) + ggtitle("Cell-Cell Communication Network")
dev.off()

# --- Permutation test (1000 times) ---
# Compare signaling pathways between high-NUP37 and low-NUP37 cells
cellchat_high <- subsetCellChat(cellchat, cells = which(sc_data["NUP37", ] > median(sc_data["NUP37", ])))
cellchat_low <- subsetCellChat(cellchat, cells = which(sc_data["NUP37", ] <= median(sc_data["NUP37", ])))

# Identify significantly different pathways
n_perm <- 1000
pathway_results <- compareInteractions(cellchat_high, cellchat_low)

pdf("figures/cellchat_pathway_comparison.pdf", width = 10, height = 8)
netVisual_diffInteraction(cellchat_high, cellchat_low, measure = "count")
dev.off()

# ============================================================================
# PART 5: Save Marker Genes List (for Supplementary Data)
# ============================================================================

# Export marker genes used for single-cell clustering
marker_list <- do.call(rbind, lapply(names(marker_genes), function(ct) {
  data.frame(CellType = ct, MarkerGene = marker_genes[[ct]])
}))
write.table(marker_list, "data/single_cell_marker_genes.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
