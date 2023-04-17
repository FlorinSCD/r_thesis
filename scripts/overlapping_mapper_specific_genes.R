# Transform the dgCMatrix m_af to a dataframe
m_af_df <- data.frame(m_af)

# Make a new Column with the sum of the rows to find out which rows sum up to 0
m_af_non_zero <-  m_af_df %>% mutate(count=rowSums(.!=0))

# Filter the dataframe to have only the rows that are non-zero
m_af_non_zero <- filter(m_af_non_zero, count > 0)

# Check which rows intesect to find out which Genes are in common
intersect(row.names(m_kb_non_zero), row.names(m_af_non_zero))

# Check which rows differ, to find out which Genes are found differently
# These values appear in m_kb but do not appear in m_af
setdiff(row.names(m_kb_non_zero_1), row.names(m_af_non_zero_1))

# By reversing the dataframes we will get different answers
# These values appear in m_af but do not appear in m_kb
setdiff(row.names(m_af_non_zero_1), row.names(m_kb_non_zero_1))

# Number of Genes per cell barcode is the number of genes expressed per barcode (DETECTED GENES / CELL)
colSums(m_af_non_zero_1 != 0)

# UMIs per cell is the number of dedulicatedReads (UMI PER CELL)
colSums(Filter(is.numeric, m_af_non_zero_1))


# Get common genes from the two dataframes and assign it to a new dataframe
common_genes <- intersect(row.names(m_kb_non_zero_1), row.names(m_af_non_zero_1))
common_genes_kb <- m_kb_non_zero_1[common_genes,]
common_genes_af <- m_af_non_zero_1[common_genes,]

common_genes_af <-common_genes_af[291:291]
common_genes_kb <-common_genes_kb[294:294]

cor(common_genes_af$count, common_genes_kb$count, method = 'spearman')
cor(common_genes_af$count, common_genes_kb$count, method="pearson")

# Plot 
plot(common_genes_af$count, common_genes_kb$count, main="Correlation between common genes af - kb", xlab="common_genes_af", ylab="common_genes_kb")

# Merged DF
merged_df <- merge(common_genes_kb, common_genes_af, by='row.names', all=TRUE)

# Plot the Merge 
ggplot(data = merged_df, mapping = aes(x=count.x, y=count.y)) 
  + geom_point(alpha = 0.1, color ='blue') 
  + geom_smooth(method='lm') 
  + labs(x='Kallisto-Bustools Counts', y='Alevin-Fry Counts', title='Alevin-Fry and Kallisto-Bustools Counts Correlation') 
  + theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))