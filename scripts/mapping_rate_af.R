
mapping_rate_af(af_file) {
  
  # Load in the file using read.table since it's a tab delimited file
  features_txt <- read.table(file(af_file), header=TRUE)
  
  # Sum of the Mapped Reads over the sum of the Corrected Reads
  x <- sum(features_txt[1:290,3:3])/sum(features_txt[1:290,2:2])
  
  return(x)
  
}
