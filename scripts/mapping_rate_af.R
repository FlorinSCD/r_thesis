# Load in the file using read.table since it's a tab delimited file
features_txt <- read.table(file("/Users/florin/Desktop/Thesis_work/pipeline_runs/20230310_sort_real_af/alevin_fry_map/SCD-TEST-s701_output/SCD-TEST-s701_feature.txt"), header=TRUE)

# Sum of the Mapped Reads over the sum of the Corrected Reads
sum(features_txt[1:290,3:3])/sum(features_txt[1:290,2:2])