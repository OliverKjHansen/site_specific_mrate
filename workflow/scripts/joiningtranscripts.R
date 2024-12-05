library("tidyverse")

# """""""
# 
# """""""


args = commandArgs(trailingOnly=TRUE)
possible_lofs <- args[1]
annotated_positions <- args[2]
output <-args[3]


df1 <-read_delim(file = possible_lofs, col_names = FALSE)
colnames(df1) <- c("chrom","pos","ref", "alt", "transcript", "type", "sequence")
df2 <-read_delim(file = annotated_positions, col_names = TRUE)

res <-inner_join(df1, df2, by = c("chrom","pos"))

write.table(res, file = output, sep='\t', quote=FALSE, row.names = FALSE)
