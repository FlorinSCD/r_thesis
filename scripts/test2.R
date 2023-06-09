library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(scico)
library(scales)

read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "CsparseMatrix")
  # The matrix read has cells in rows
  genes <- readLines(file(paste0(dir, "/", "features.tsv")))
  barcodes <- readLines(file(paste0(dir, "/", "barcodes.tsv")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

res_mat <- read_count_output("/Users/florin/Desktop/Thesis_work/Matrices/Sort/Real/kb-s701/SCD-TEST-s701-filtered-feature-bc-matrix", "matrix")
dim(res_mat)

seu <- CreateSeuratObject(res_mat, min.cells = 3, min.features = 200)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# NORMALIZING THE DATA 

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu <- NormalizeData(seu)

# IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# SCALING THE DATA 

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

# PERFORM LINEAR DIMENSIONAL REDUCTION

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

# Examine and visualize PCA results a few different ways
print(seu[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seu, dims = 1:2, reduction = "pca")

DimPlot(seu, reduction = "pca")

DimHeatmap(seu, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE)

# DETERMINE THE 'DIMENSIONALITY' OF THE DATASET

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
seu <- JackStraw(seu, num.replicate = 100)
seu <- ScoreJackStraw(seu, dims = 1:20)

JackStrawPlot(seu, dims = 1:15)

ElbowPlot(seu)

# CLUSTER THE CELLS

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seu), 5)

# RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
seu <- RunUMAP(seu, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seu, reduction = "umap")

saveRDS(seu, file = "../output/pbmc_tutorial.rds")

# FINDING DIFFERENTIALLY EXPRESSED FEATURES (CLUSTER BIOMARKERS)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(seu, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 2s from clusters 0 and 1
cluster5.markers <- FindMarkers(seu, ident.1 = 2, ident.2 = c(0, 1), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seu.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(seu, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(seu, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(seu, features = c("NKG7", "LEF1"), slot = "counts", log = TRUE)


FeaturePlot(seu, features = c("MS4A1", "GNLY", "CD3E", "FCGR3A", "LYZ", "CD8A"))


seu.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seu, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T")
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(seu, file = "../output/pbmc3k_final.rds")
