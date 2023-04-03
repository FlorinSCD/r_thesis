library(Matrix)
library(stringr)

# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name, kb) {
  
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "CsparseMatrix")
  
  # The matrix read has cells in rows
  genes_df <- readLines(file(paste0(dir, "/", "features.tsv")))
  barcodes <- readLines(file(paste0(dir, "/", "barcodes.tsv")))
  
  df <- read.csv(file("/Users/florin/Desktop/Thesis_work/r_project/data/gene_mapping_modified.csv"))
  
  if(kb) {
    # Remove the version number of the ensembl gene IDs
    # This is only needed for Kallisto-Bustools. Alevin-Fry does not produce the output in a way that this is needed
    genes_df <- str_replace(genes_df, pattern = ".[0-9]+$", replacement = "")
  }
  
  # Make genes vector into a data frame, rename genes column and create new column
  genes_df <- data.frame(genes_df)
  colnames(genes_df)[1] <- 'gene_ids'
  genes_df['genes'] <- ''
  
  # This is to replace "\" in the csv file
  for(i in 1:nrow(df)) {
    for(j in 1:ncol(df)) {
      vec <- gsub("\"", "", df[i, j])
      df[i,j] <- vec
    }
  }
  
  # To merge csv_df genes column onto genes_df genes column
  genes_df$genes <- df$genes[match(genes_df$gene_ids, df$gene_ids)]
  
  # Make genes column from data frame into a vector 
  genes_vector <- genes_df$genes
  
  colnames(m) <- barcodes
  rownames(m) <- genes_vector
  return(m)
}