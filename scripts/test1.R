library(Matrix)
library(Seurat)
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
  return(genes)
  
}

res_mat <- read_count_output("/Users/florin/Desktop/Thesis_work/Matrices/Sort/Real/kb-s701/SCD-TEST-s701-filtered-feature-bc-matrix", "matrix")
dim(res_mat)

seu <- CreateSeuratObject(res_mat, min.cells = 3, min.features = 200)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
# Visualize QC metrics as a violin plot
options(repr.plot.width=12, repr.plot.height=6)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

options(repr.plot.width=9, repr.plot.height=6)
ggplot(seu@meta.data, aes(nCount_RNA, nFeature_RNA)) +
  geom_hex(bins = 100) +
  scale_fill_scico(palette = "devon", direction = -1, end = 0.9) +
  scale_x_log10(breaks = breaks_log(12)) + 
  scale_y_log10(breaks = breaks_log(12)) + annotation_logticks() +
  labs(x = "Total UMI counts", y = "Number of genes detected") +
  theme(panel.grid.minor = element_blank())

seu <- subset(seu, subset = percent.mt < 3)

# seu <- NormalizeData(seu) %>% ScaleData()

# ANALYSIS

# Identify highly variable genes

seu <- FindVariableFeatures(seu, nfeatures = 3000)
top10 <- head(VariableFeatures(seu), 10)
plot1 <- VariableFeaturePlot(seu, log = FALSE)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# PCA

seu <- RunPCA(seu, verbose = FALSE, npcs = 20) # uses HVG by default
ElbowPlot(seu, ndims = 20)

PCAPlot(seu)

# Clustering and visualization

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu)
PCAPlot(seu)


# tSNE

seu <- RunTSNE(seu, dims = 1:10)
TSNEPlot(seu)

# UMAP

seu <- RunUMAP(seu, dims = 1:10, verbose = FALSE)
UMAPPlot(seu)

