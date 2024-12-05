library("tidyverse")

args = commandArgs(trailingOnly=TRUE)
genovo <- args[1]
canconical <- args[2]
wtranscripts <- args[3]
haplo <- args[4]
output <- args[5]

#observed <- read.delim(file = genovo) %>% filter(mutation_type == "LoF" | mutation_type == "splice_site"| mutation_type == "nonsense"| mutation_type =="frameshift_indel") %>% group_by(region) %>% summarise(obs = sum(observed))
observed <- read.delim(file = genovo) %>% filter(mutation_type == "LoF") %>% group_by(region) %>% summarise(obs = sum(observed), exp = sum(expected), exp_up = sum(expected_upper))
expected <- read.delim(file = wtranscripts) %>% group_by(transcript) %>% summarise(se = sum(prob_1se), min = sum(prob_min))

canonical <- read.delim(file = canconical) %>% filter(Ensembl.Canonical == "1")
#canonical$Transcript.stable.ID.version
#canonical$NCBI.gene..formerly.Entrezgene..accession
#canonical$Gene.name
#canonical$Gene.stable.ID

res <- expected %>% full_join(observed, by = c("transcript"="region")) %>% filter(!is.na(obs) & !is.na(min))
#res2 <- observed %>% full_join(expected, by = c("region"="transcript")) %>% filter(!is.na(obs) & !is.na(se))
non_matched <- expected %>% full_join(observed, by = c("transcript"="region")) %>% filter(is.na(obs) | is.na(se))
#non_matched1 <- observed %>% full_join(expected, by = c("region"="transcript")) %>% filter(is.na(obs) | is.na(se))
print(non_matched)

final <- canonical %>% full_join(res, by = c("Transcript.stable.ID.version"="transcript")) %>% filter(!is.na(obs) & !is.na(Gene.name)) %>% mutate(min_OE=obs/min,
                                                                                                                                                   se_OE=obs/se)
non_matched_final <- canonical %>% full_join(res, by = c("Transcript.stable.ID.version"="transcript")) %>% filter(is.na(obs) | is.na(Gene.name))
print(non_matched_final)

curated <- read_delim(file = haplo, delim = "\t") %>% filter(`Haploinsufficiency Score` == 3) %>% separate(`Genomic Location`,sep = ":",into = c("chr","pos")) 

final$haploinsufficiency <- final$NCBI.gene..formerly.Entrezgene..ID %in% curated$`Gene ID`

##add loeuf analysis

write_delim(final, file = output, col_names = TRUE, delim = "\t")