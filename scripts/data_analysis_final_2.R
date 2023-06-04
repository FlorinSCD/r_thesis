library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpointdensity)
library(scico)
library(scales)
library(cowplot)

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
  feature_plot_all <- feature_plot_1 + feature_plot_2 + feature_plot_3
  
  feature_plot_1_ggplot <- ggplot(seu[[]], aes(percent.mt, nCount_RNA)) + geom_point() + scale_y_log10()
  feature_plot_2_ggplot <- ggplot(seu[[]], aes(percent.rp, nCount_RNA)) + geom_point() + scale_y_log10()
  feature_plot_3_ggplot <- ggplot(seu[[]], aes(nFeature_RNA, nCount_RNA)) + geom_point() + scale_y_log10()
  feature_plot_all_ggplot <- feature_plot_1_ggplot + feature_plot_2_ggplot + feature_plot_3_ggplot
  
  # UMIs per cell
  umis_per_cell <- ggplot(seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of UMIs per Cell of the g001 Dataset") + ylab("Count") + theme(plot.title = element_text(hjust = 0.5))
  
  # Genes per cell
  genes_per_cell <- ggplot(seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of Genes per Cell of the g001 Dataset") + ylab("Count") + theme(plot.title = element_text(hjust = 0.5))
  
  # Mitochondrial genes per cell
  mitochondrial_genes_per_cell <- ggplot(seu[[]], aes (percent.mt)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Mitochondrial Genes per Cell of the g001 Dataset") + xlab("Mitochondrial Genes Percentage") + ylab("Count") + theme(plot.title = element_text(hjust = 0.5))
  
  # Ribosomal genes per cell
  ribosomal_genes_per_cell <- ggplot(seu[[]], aes (percent.rp)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Ribosomal Genes per Cell of the g001 Dataset") + xlab("Ribosomal Genes Percentage") + ylab("Count") + theme(plot.title = element_text(hjust = 0.5))
  
  # TOTAL UMI COUNTS PER PERCENTAGE MITOCHONDRIAL
  umi_vs_mito <- ggplot(seu@meta.data, aes(percent.mt, nCount_RNA)) +
    geom_pointdensity() +
    scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
    labs(x = "Percentage Mitochondrial", y = "Total UMI Counts") +
    ggtitle("Total UMI Counts per Percentage Mitochondrial of g001 Dataset") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_log10()
  
  # PROCESSING
  
  # Normalize the data
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Scale the data
  all.genes <- rownames(seu)
  seu <- ScaleData(seu, features = all.genes)
  
  # IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seu), 10)
  plot1 <- VariableFeaturePlot(seu, log = FALSE)
  LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  # Scaling the data
  all.genes <- rownames(seu)
  seu <- ScaleData(seu, features = all.genes)
  
  # Perform Linear Dimensional reduction
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  elbow_plot <- ElbowPlot(seu, ndims = 20)
  
  print(seu[["pca"]], dims = 1:5, nfeatures = 5)
  
  # Clustering and visualization
  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu)
  pca_plot <- PCAPlot(seu)
  
  # tSNE
  seu <- RunTSNE(seu, dims = 1:10)
  tsne_plot <- TSNEPlot(seu) + ggtitle("tSNE Plot of the g001 Dataset") + theme(plot.title = element_text(hjust = 0.5))
  
  # UMAP
  seu <- RunUMAP(seu, dims = 1:10, verbose = FALSE)
  umap_plot <- UMAPPlot(seu) + ggtitle("UMAP Plot of the g001 Dataset") + theme(plot.title = element_text(hjust = 0.5))
  #umap_plot <- DimPlot(object = seu_combined, reduction = 'tsne', group.by = 'orig.ident', cols = c('ss1' = 'indianred1', 'ss2' = 'royalblue1', 'ss3' = 'lightgreen')) + ggtitle("tSNE of all Datasets Combined") + theme(plot.title = element_text(hjust = 0.5))
  
  
  # CLUSTERING
  
  #clustered_umap <- DimPlot(seu, reduction = "umap")  + ggtitle("Clustered UMAP") + theme(plot.title = element_text(hjust = 0.5))
  #clustered_tsne <- DimPlot(seu, reduction = "tsne")  + ggtitle("Clustered tSNE") + theme(plot.title = element_text(hjust = 0.5))
  
  # DE ON CLUSERS
  
  seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  seu.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
  # Heatmap of top 10 genes per cluster
  
  top10 <- seu.markers %>% group_by(cluster) %>% top_n(n = 10)
  heat_map_all <- DoHeatmap(seu, features = top10$gene) + ggtitle("Heatmap of Top 10 Genes per Cluster") + theme(plot.title = element_text(hjust = 0.5)) #  + NoLegend() 
  
  top_gene_clusters <- list()
  
  for (i in 1:(length(top10)-1)) {
    top_gene_clusters[[i]] <- DotPlot(seu, features = top10$gene[(1+10*(i-1)):(10+10*(i-1))]) + RotatedAxis() + coord_flip() + ggtitle(paste("Top 10 Genes For Cluster ", i))
  }
  
  top_gene_clusters_all <- plot_grid(plotlist = top_gene_clusters) #+ ggtitle("Average Expression and Percentage Expressed of Top Genes for each Cluster")
  
  return(list(violin_plot_1, feature_plot_all, feature_plot_all_ggplot, umis_per_cell, genes_per_cell, mitochondrial_genes_per_cell, ribosomal_genes_per_cell, umi_vs_mito, elbow_plot, pca_plot, tsne_plot, umap_plot, clustered_umap, clustered_tsne, heat_map_all, top_gene_clusters_all))
  #return(list(violin_plot_1, feature_plot_1, feature_plot_2, feature_plot_3, feature_plot_1_ggplot, umis_per_cell, genes_per_cell,  mitochondrial_genes_per_cell, ribosomal_genes_per_cell, umi_vs_mito, elbow_plot, pca_plot, tsne_plot, umap_plot))
  
}

combine_datasets <- function(m1, m2, m3, m4) {
  
  m1_seu <- CreateSeuratObject(counts = m1, project = 'ss1')
  m2_seu <- CreateSeuratObject(counts = m2, project = 'ss2')
  m3_seu <- CreateSeuratObject(counts = m3, project = 'ss3')
  m4_seu <- CreateSeuratObject(counts = m4, project = 'ss4')
  
  
  seu_combined <- merge(m1_seu, y = c(m2_seu, m3_seu, m4_seu), add.cell.ids = c("ss1", "ss2", "ss3", "ss4"))
  
  m1_seu[['percent.mt']] <- PercentageFeatureSet(m1_seu, pattern = "^MT")
  m1_seu[["percent.rp"]] <- PercentageFeatureSet(m1_seu, pattern = "^RP")
  
  m2_seu[['percent.mt']] <- PercentageFeatureSet(m2_seu, pattern = "^MT")
  m2_seu[["percent.rp"]] <- PercentageFeatureSet(m2_seu, pattern = "^RP")
  
  m3_seu[['percent.mt']] <- PercentageFeatureSet(m3_seu, pattern = "^MT")
  m3_seu[["percent.rp"]] <- PercentageFeatureSet(m3_seu, pattern = "^RP")
  
  m4_seu[['percent.mt']] <- PercentageFeatureSet(m4_seu, pattern = "^MT")
  m4_seu[["percent.rp"]] <- PercentageFeatureSet(m4_seu, pattern = "^RP")
  
  # QUALITY METIRCS
  
  # UMIs per cell
  
  umis_per_cell <- ggplot(m4_seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("UMI Count per Cell for Cell Ranger") + xlab("Number of Molecules Detected") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umis_per_cell1 <- ggplot(m1_seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("UMI Count per Cell for Alevin-Fry") + ylab("Count") + xlab("Number of Molecules Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umis_per_cell2 <- ggplot(m2_seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("UMI Count per Cell for Kallisto Bustools") + ylab("Count") + xlab("Number of Molecules Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umis_per_cell3 <- ggplot(m3_seu[[]], aes (nCount_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("UMI Count per Cell for STARsolo") + ylab("Count") + xlab("Number of Molecules Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umis_per_cell_combined <- umis_per_cell1 + umis_per_cell2 + umis_per_cell3 + umis_per_cell
  ggsave("results/umis_per_cell_combined.png", umis_per_cell_combined, width = 2560, height = 1369, units = c("px"))
  
  # Genes per cell
  genes_per_cell <- ggplot(m4_seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of Genes per Cell for Cell Ranger") + xlab("Number of Genes Detected") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  genes_per_cell1 <- ggplot(m1_seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of Genes per Cell for Alevin-Fry") + ylab("Count") + xlab("Number of Genes Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  genes_per_cell2 <- ggplot(m2_seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of Genes per Cell for Kallisto Bustools") + ylab("Count") + xlab("Number of Genes Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  genes_per_cell3 <- ggplot(m3_seu[[]], aes (nFeature_RNA)) + geom_histogram(bins = 80) + scale_x_log10() + ggtitle("Count of Genes per Cell for STARsolo") + ylab("Count") + xlab("Number of Genes Detected") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  genes_per_cell_combined <- genes_per_cell1 + genes_per_cell2 + genes_per_cell3 + genes_per_cell
  ggsave("results/genes_per_cell_combined.png", genes_per_cell_combined, width = 2560, height = 1369, units = c("px"))
  
  # Mitochondrial genes per cell
  mitochondrial_genes_per_cell <- ggplot(m4_seu[[]], aes (percent.mt)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Mitochondrial Genes per Cell for Cell Ranger") + xlab("Percentage of Mitochondrial Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  mitochondrial_genes_per_cell1 <- ggplot(m1_seu[[]], aes (percent.mt)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Mitochondrial Genes per Cell for Alevin-Fry") + xlab("Percentage of Mitochondrial Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  mitochondrial_genes_per_cell2 <- ggplot(m2_seu[[]], aes (percent.mt)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Mitochondrial Genes per Cell for Kallisto Bustools") + xlab("Percentage of Mitochondrial Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  mitochondrial_genes_per_cell3 <- ggplot(m3_seu[[]], aes (percent.mt)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Mitochondrial Genes per Cell for STARsolo") + xlab("Percentage of Mitochondrial Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  mitochondrial_genes_per_cell_combined <- mitochondrial_genes_per_cell1 + mitochondrial_genes_per_cell2 + mitochondrial_genes_per_cell3 + mitochondrial_genes_per_cell 
  ggsave("results/mitochondrial_genes_per_cell_combined.png", mitochondrial_genes_per_cell_combined, width = 2560, height = 1369, units = c("px"))
  
  # Ribosomal genes per cell
  ribosomal_genes_per_cell <- ggplot(m4_seu[[]], aes (percent.rp)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Ribosomal Genes per Cell for Cell Ranger") + xlab("Percentage of Ribosomal Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ribosomal_genes_per_cell1 <- ggplot(m1_seu[[]], aes (percent.rp)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Ribosomal Genes per Cell For Alevin-Fry") + xlab("Percentage of Ribosomal Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ribosomal_genes_per_cell2 <- ggplot(m2_seu[[]], aes (percent.rp)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Ribosomal Genes per Cell for Kallisto Bustools") + xlab("Percentage of Ribosomal Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ribosomal_genes_per_cell3 <- ggplot(m3_seu[[]], aes (percent.rp)) + geom_histogram(bins = 80) + xlim(0, 100) + ggtitle("Count of Ribosomal Genes per Cell for STARsolo") + xlab("Percentage of Ribosomal Genes") + ylab("Count") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ribosomal_genes_per_cell_combined <- ribosomal_genes_per_cell1 + ribosomal_genes_per_cell2 + ribosomal_genes_per_cell3 + ribosomal_genes_per_cell
  ggsave("results/ribosomal_genes_per_cell_combined.png", ribosomal_genes_per_cell_combined, width = 2560, height = 1369, units = c("px"))
  
  # TOTAL UMI COUNTS PER PERCENTAGE MITOCHONDRIAL
  umi_vs_mito1 <- ggplot(m4_seu@meta.data, aes(percent.mt, nCount_RNA)) +
    geom_pointdensity() +
    scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
    labs(x = "Percentage of Mitochondrial Cells", y = "Total UMI Counts") +
    ggtitle("Total UMI Counts per Percentage Mitochondrial\nsfor Cell Ranger") +
    theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
    scale_y_log10()
  umi_vs_mito2 <- ggplot(m1_seu@meta.data, aes(percent.mt, nCount_RNA)) +
    geom_pointdensity() +
    scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
    labs(x = "Percentage of Mitochondrial Cells", y = "Total UMI Counts") +
    ggtitle("Total UMI Counts per Percentage Mitochondrial\nfor Alevin-Fry") +
    theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
    scale_y_log10()
  umi_vs_mito3 <- ggplot(m2_seu@meta.data, aes(percent.mt, nCount_RNA)) +
    geom_pointdensity() +
    scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
    labs(x = "Percentage of Mitochondrial Cells", y = "Total UMI Counts") +
    ggtitle("Total UMI Counts per Percentage Mitochondrial\nfor Kallisto Bustools") +
    theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
    scale_y_log10()
  umi_vs_mito <- ggplot(m3_seu@meta.data, aes(percent.mt, nCount_RNA)) +
    geom_pointdensity() +
    scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
    labs(x = "Percentage of Mitochondrial Cells", y = "Total UMI Counts") +
    ggtitle("Total UMI Counts per Percentage Mitochondrial\nfor STARsolo") +
    theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
    scale_y_log10()
  umi_vs_mito_combined <- umi_vs_mito1 + umi_vs_mito2 + umi_vs_mito3 + umi_vs_mito
  ggsave("results/umi_vs_mito_combined.png", umi_vs_mito_combined, width = 2560, height = 1369, units = c("px"))
  
  # PROCESSING
  
  # Normalize the data
  seu_combined <- NormalizeData(seu_combined, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Scale the data
  all.genes <- rownames(seu_combined)
  seu_combined <- ScaleData(seu_combined, features = all.genes)
  
  # IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION)
  seu_combined <- FindVariableFeatures(seu_combined, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seu_combined), 10)
  plot1 <- VariableFeaturePlot(seu_combined, log = FALSE)
  LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  # Scaling the data
  all.genes <- rownames(seu_combined)
  seu_combined <- ScaleData(seu_combined, features = all.genes)
  
  # Perform Linear Dimensional reduction
  seu_combined <- RunPCA(seu_combined, features = VariableFeatures(object = seu_combined))
  elbow_plot <- ElbowPlot(seu_combined, ndims = 20) + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ggsave("results/elbow_plot.jpeg", elbow_plot, width = 2560, height = 1369, units = c("px"))
  
  print(seu_combined[["pca"]], dims = 1:5, nfeatures = 5)
  
  # Clustering and visualization
  seu_combined <- FindNeighbors(seu_combined, dims = 1:10)
  seu_combined <- FindClusters(seu_combined)
  pca_plot <- PCAPlot(seu_combined) + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) 
  ggsave("results/pca_plot.jpeg", pca_plot, width = 2560, height = 1369, units = c("px"))
  
  # tSNE per Library
  seu_combined <- RunTSNE(seu_combined, dims = 1:10)
  tsne_plot <- DimPlot(object = seu_combined, reduction = 'tsne', group.by = 'orig.ident', cols = c('ss1' = 'indianred1', 'ss2' = 'grey87', 'ss3' = 'grey87', 'ss4' = 'grey87')) + ggtitle("tSNE Plot for Cell Ranger") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  tsne1_plot <- DimPlot(object = seu_combined, reduction = 'tsne', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = adjustcolor('royalblue1', alpha.f = 1), 'ss3' = adjustcolor('grey87', alpha.f = 1), 'ss4' = 'grey87')) + ggtitle("tSNE Plot for Alevin-Fry") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  tsne2_plot <- DimPlot(object = seu_combined, reduction = 'tsne', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = 'grey87', 'ss3' = 'lightgreen', 'ss4' = 'grey87')) + ggtitle("tSNE Plot for Kallisto Bustools") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  tsne3_plot <- DimPlot(object = seu_combined, reduction = 'tsne', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = 'grey87', 'ss3' = 'grey87','ss4' = 'purple')) + ggtitle("tSNE Plot for STARsolo") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  tsne_all <- tsne1_plot + tsne2_plot + tsne3_plot + tsne_plot
  ggsave("results/tsne_all.png", tsne_all, width = 2560, height = 1369, units = c("px"))
  
  # UMAP per Library
  seu_combined <- RunUMAP(seu_combined, dims = 1:10, verbose = FALSE)
  umap_plot <- DimPlot(object = seu_combined, reduction = 'umap', group.by = 'orig.ident', cols = c('ss1' = 'indianred1', 'ss2' = 'grey87', 'ss3' = 'grey87', 'ss4' = 'grey87')) + ggtitle("UMAP Plot for Cell Ranger") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umap_ss1 <- DimPlot(object = seu_combined, reduction = 'umap', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = adjustcolor('royalblue1', alpha.f = 1), 'ss3' = adjustcolor('grey87', alpha.f = 1), 'ss4' = 'grey87')) + ggtitle("UMAP Plot for Alevin-Fry") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umap_ss2 <- DimPlot(object = seu_combined, reduction = 'umap', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = 'grey87', 'ss3' = 'lightgreen', 'ss4' = 'grey87')) + ggtitle("UMAP Plot for Kallisto Bustools") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umap_ss3 <- DimPlot(object = seu_combined, reduction = 'umap', group.by = 'orig.ident', cols = c('ss1' = 'grey87', 'ss2' = 'grey87', 'ss3' = 'grey87','ss4' = 'purple'))  + ggtitle("UMAP Plot for STARsolo") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  umap_all <- umap_ss1 + umap_ss2 + umap_ss3 + umap_plot
  ggsave("results/umap_all.png", umap_all, width = 2560, height = 1369, units = c("px"))
  
  # CLUSTERING
  
  clustered_umap <- DimPlot(seu_combined, reduction = "umap")  + ggtitle("Clustered UMAP") + xlab("UMAP 1") + ylab("UMAP 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  clustered_tsne <- DimPlot(seu_combined, reduction = "tsne")  + ggtitle("Clustered tSNE") + xlab("tSNE 1") + ylab("tSNE 2") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  ggsave("results/clustered_umap.png", clustered_umap, width = 2560, height = 1369, units = c("px"))
  ggsave("results/clustered_tsne.png", clustered_tsne, width = 2560, height = 1369, units = c("px"))
  
  # DE ON CLUSTERS
  
  # find all markers of cluster 2
  #cluster2.markers <- FindMarkers(seu_combined, ident.1 = 2, min.pct = 0.25)
  #head_2 <- head(cluster2.markers, n = 5)
  
  #node <- list()
  
  #for (i in 1:10) {
  
  #  cluster_marker <- FindMarkers(seu_combined, ident.1 = i, min.pct = 0.25)
  #  node[[i]] <- head(cluster_marker, n = 7)
  
  #}
  
  #DefaultAssay(seu_combined) <- "RNA"
  #nk.markers <- FindConservedMarkers(seu_combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
  #head(nk.markers)
  
  seu_combined.markers <- FindAllMarkers(seu_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  seu_combined.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
  # Heatmap of top 10 genes per cluster
  
  top10 <- seu_combined.markers %>% group_by(cluster) %>% top_n(n = 10)
  heat_map_all <- DoHeatmap(seu_combined, features = top10$gene) + ggtitle("Heatmap of Top 10 Genes per Cluster") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 5)) #  + NoLegend() 
  ggsave("results/heat_map_all.png", heat_map_all, width = 2560, height = 1369, units = c("px"))
  
  top_gene_clusters <- list()
  
  for (i in 1:(length(top10)-1)) {
    top_gene_clusters[[i]] <- DotPlot(seu_combined, features = top10$gene[(1+10*(i-1)):(10+10*(i-1))]) + RotatedAxis() + coord_flip() + theme(text = element_text(size = 8), plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 5), legend.key.height = unit(0.15, 'cm'), legend.key.width = unit(0.05, 'cm')) #+ NoLegend() 
  }
  
  top_gene_clusters_all <- plot_grid(plotlist = top_gene_clusters) + ggtitle("Average Expression and Percentage Expressed of Top Genes for each Cluster") + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) #  + NoLegend() 
  ggsave("results/top_gene_clusters_all.jpeg", top_gene_clusters_all, width = 3000, height = 1369, units = c("px"))
  
  #return(list(umis_per_cell,genes_per_cell,mitochondrial_genes_per_cell,ribosomal_genes_per_cell, umi_vs_mito, elbow_plot, pca_plot, tsne_plot, umap_plot, umap_plot_2))
  return(list(umis_per_cell_combined, genes_per_cell_combined, mitochondrial_genes_per_cell_combined, ribosomal_genes_per_cell_combined, umi_vs_mito_combined, elbow_plot, pca_plot, tsne_all, umap_all, clustered_umap, clustered_tsne, heat_map_all, top_gene_clusters_all))
  
}

# Mapping Rate of Kallisto-Bustools

mapping_rate_kb<- function(kb_file) {
  
  # run_info.json
  
  # Load in the file using fromJSON since it's a json file
  library("rjson")
  
  # run_info.json
  features_kb <- fromJSON(file=kb_file)
  
  # Load the p_pseudoaligned 
  features_kb[6]
  
  # In case needed to know, this creates a dataframe of the json file
  features_kb <- as.data.frame(features_kb)
  
  # Load the actual value of p_pseudoaligned
  p_pseudoaligned <- features_kb$p_pseudoaligned
  
  return(p_pseudoaligned)
  
}

# Mapping Rate of Alevin-Fry

mapping_rate_af<- function(af_file) {
  
  # SCD-WP-g001_feature.txt
  
  # Load in the file using read.table since it's a tab delimited file
  features_txt <- read.table(file(af_file), header=TRUE)
  
  # Sum of the Mapped Reads over the sum of the Corrected Reads
  alevin_fry_mapping_rate <- sum(features_txt[1:290,3:3])/sum(features_txt[1:290,2:2])
  
  return(alevin_fry_mapping_rate)
  
}

