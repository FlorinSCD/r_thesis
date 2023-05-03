
mapping_rate_kb(kb_file) {

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