library(Matrix)
library(stringr)

load_truth_mtx <- function(dir) {
  
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/trueUMIcounts.mtx"))
  m <- Matrix::t(m)
  m <- as(m, "CsparseMatrix")
  
  return(m)
  
}
