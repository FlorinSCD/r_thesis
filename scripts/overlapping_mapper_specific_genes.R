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

# Number of Genes per cell barcode is the number of genes expressed per barcode
colSums(m_af_non_zero_1 != 0)

# UMIs per cell is the number of dedulicatedReads
colSums(Filter(is.numeric, m_af_non_zero_1))
