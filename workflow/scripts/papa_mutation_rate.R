library("tidyverse")

args = commandArgs(trailingOnly=TRUE)
input <- args[1]
snv_out <- args[2]
indel_out <- args[3]

file  <- read_delim(file = input, col_names = FALSE)
colnames(file) <- c("type", "pattern", "rate")

selected_types <- c("A2C", "A2T", "A2G", "C2G", "C2A", "C2T")
snv <- file %>% filter(type %in% selected_types) %>% mutate(type = gsub("2", "->", type))

indel <- file %>%
  filter(type %in% c("insertion", "deletion")) %>%
  pivot_wider(
    names_from = type,
    values_from = rate,
    values_fill = 0
  ) %>% mutate( frameshift = insertion+deletion, #dobbelttjek this again
                inframe = 0) %>% select(c(pattern, frameshift, inframe))

write_delim(snv, file = snv_out, col_names = FALSE, delim = " ")
write_delim(indel, file = indel_out, col_names = FALSE, delim = " ")