library("tidyverse")

args = commandArgs(trailingOnly=TRUE)
inputfile <- args[1]
output_file_1se <- args[2]
output_file_min <- args[3]

df <- read_delim(file =inputfile)
hd <- c("chrom","pos","prob_min","prob_1se","transcript_id", "mutationtype")
colnames(df) <- hd
transcripts_min <- df %>% 
  group_by(transcript_id) %>% 
  summarise(lambda_min = sum(prob_min))

transcripts_1se <- df %>% 
  group_by(transcript_id) %>% 
  summarise(lambda_1se = sum(prob_1se))

write_delim(transcripts_min, file = output_file_min, col_names = TRUE, delim = "\t")
write_delim(transcripts_1se, file = output_file_1se, col_names = TRUE, delim = "\t")