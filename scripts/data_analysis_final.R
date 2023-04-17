library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpointdensity)
library(scico)
library(scales)

data_analysis <- function(sparse_matrix) {
  
  m <- sparse_matrix
  
  # Create a Seurat Object
  seu <- CreateSeuratObject(m)
  
  # Add mitochondrial and ribosomal columns to Seurat Object
  seu[['percent.mt']] <- PercentageFeatureSet(seu, pattern = "^MT")
  seu[["percent.rp"]] <- PercentageFeatureSet(seu, pattern = "^RP")
  
  # Make a violin plot of the mitochondrial, ribosomal and nCount/nFeature_RNA data
  violin_plot_1 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  feature_plot_1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
  feature_plot_2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.rp")
  feature_plot_3 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  feature_plot_1 + feature_plot_2 + feature_plot_3
  
  feature_plot_1_ggplot <- ggplot(seu[[]], aes (percent.mt, nCount_RNA)) + geom_point() + scale_y_log10()
  
  
  
  # UMIs per cell
  umis_per_cell <- ggplot(seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 30) + scale_x_log10()
  
  # Genes per cell
  genes_per_cell <- ggplot(seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 30) + scale_x_log10()
  
  # Mitochondrial genes per cell
  mitochondrial_genes_per_cell <- ggplot(seu[[]], aes (percent.mt)) + geom_histogram(bins = 30) + xlim(0, 100)
  
  # Ribosomal genes per cell
  ribosomal_genes_per_cell <- ggplot(seu[[]], aes (percent.rp)) + geom_histogram(bins = 30) + xlim(0, 100)
  
  # TOTAL UMI COUNTS PER PERCENTAGE MITOCHONDRIAL
  umi_vs_mito <- ggplot(seu@meta.data, aes(percent.mt, nCount_RNA)) +
    geom_pointdensity() +
    scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
    labs(x = "Percentage Mitochondrial", y = "Total UMI Counts") +
    scale_y_log10()
  
  
  
  # Normalize the data
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Scale the data
  all.genes <- rownames(seu)
  seu <- ScaleData(seu, features = all.genes)
  
  # IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seu), 10)
  
  # Scaling the data
  all.genes <- rownames(seu)
  seu <- ScaleData(seu, features = all.genes)
  
  # Perform Linear Dimensional reduction
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  
  print(seu[[pca]], dims = 1:5, nfeatures = 5)
  
  
  return(list(violin_plot_1, feature_plot_1, feature_plot_2, feature_plot_3, feature_plot_1_ggplot, umis_per_cell, genes_per_cell,  mitochondrial_genes_per_cell, ribosomal_genes_per_cell, umi_vs_mito))
  
}