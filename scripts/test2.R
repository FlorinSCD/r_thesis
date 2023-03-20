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
