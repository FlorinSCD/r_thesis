library(Matrix)
library(stringr)

# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here

load_star <- function(dir) {
  
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/matrix.mtx"))
  m <- Matrix::t(m)
  m <- as(m, "CsparseMatrix")
  
  # The matrix read has cells in rows
  genes_df <- read.table(file(paste0(dir, "/", "features.tsv")))
  barcodes <- read.table(file(paste0(dir, "/", "barcodes.tsv")))
  
  genes_vector <- genes_df$V2
  barcodes <- barcodes$V1
  
  m <- t(m)
  
  colnames(m) <- barcodes
  rownames(m) <- genes_vector
  
  return(m)
  
}

load_af <- function(dir) {
  
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/matrix.mtx"))
  m <- Matrix::t(m)
  m <- as(m, "CsparseMatrix")
  
  # The matrix read has cells in rows
  genes_df <- read.table(file(paste0(dir, "/", "features.tsv")))
  barcodes <- read.table(file(paste0(dir, "/", "barcodes.tsv")))
  
  gene_mapping <- read.csv(file("/Users/florin/Desktop/Thesis_work/r_project/data/gene_mapping_modified.csv"))
  
  # Make genes vector into a data frame, rename genes column and create new column
  genes_df <- data.frame(genes_df)
  colnames(genes_df)[1] <- 'gene_ids'
  genes_df['genes'] <- ''
  
  m_fun <- function(x) {
    
    return (x <- gsub("\"", "", x))
    
  }
  
  gene_mapping <- data.frame(lapply(gene_mapping, m_fun))
  
  # To merge csv_df genes column onto genes_df genes column
  genes_df$genes <- gene_mapping$genes[match(genes_df$gene_ids, gene_mapping$gene_ids)]
  
  # Make genes column from data frame into a vector 
  genes_vector <- genes_df$genes
  
  barcodes <- barcodes$V1
  
  colnames(m) <- barcodes
  rownames(m) <- genes_vector
  
  return(m)
  
}

load_kb <- function(dir) {
  
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/matrix.mtx"))
  m <- Matrix::t(m)
  m <- as(m, "CsparseMatrix")
  
  # The matrix read has cells in rows
  genes_df <- read.table(file(paste0(dir, "/", "features.tsv")))
  barcodes <- read.table(file(paste0(dir, "/", "barcodes.tsv")))
  
  gene_mapping <- read.csv(file("/Users/florin/Desktop/Thesis_work/r_project/data/gene_mapping_modified.csv"))
  
  # Remove the version number of the ensembl gene IDs
  # This is only needed for Kallisto-Bustools. Alevin-Fry does not produce the output in a way that this is needed
  genes_df <- str_replace(genes_df$V1, pattern = ".[0-9]+$", replacement = "")
  
  # Make genes vector into a data frame, rename genes column and create new column
  genes_df <- data.frame(genes_df)
  colnames(genes_df)[1] <- 'gene_ids'
  genes_df['genes'] <- ''
  
  m_fun <- function(x) {
    
    return (x <- gsub("\"", "", x))
    
  }
  
  gene_mapping <- data.frame(lapply(gene_mapping, m_fun))
  
  # To merge csv_df genes column onto genes_df genes column
  genes_df$genes <- gene_mapping$genes[match(genes_df$gene_ids, gene_mapping$gene_ids)]
  
  # Make genes column from data frame into a vector 
  genes_vector <- genes_df$genes
  
  barcodes <- barcodes$V1
  
  colnames(m) <- barcodes
  rownames(m) <- genes_vector
  
  return(m)
  
}