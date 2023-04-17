# Load in the file using fromJSON since it's a json file
library("rjson")
features_kb <- fromJSON(file="/Users/florin/Desktop/Thesis_work/pipeline_runs/20230310_sort_real_kb/kallisto_bustools_map/SCD-TEST-s703_output/run_info.json")

# Load the p_pseudoaligned 
features_kb[6]

# In case needed to know, this creates a dataframe of the json file
features_kb <- as.data.frame(features_kb)

# Load the actual value of p_pseudoaligned
features_kb$p_pseudoaligned