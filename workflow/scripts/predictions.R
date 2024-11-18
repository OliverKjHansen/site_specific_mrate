library(tidyverse)
#library(stringi)
library(glmnet)

# inputfile
args = commandArgs(trailingOnly=TRUE)

file <- args[1]
model <- args[2]
levelsval <- args[3]
log_model <-args[4]
output <- args[5]

# loading in prediction dataset
df <- read.table(file, header = TRUE)
df[c("meth1")][is.na(df[c("meth1")])] <- 0
df[c("meth2")][is.na(df[c("meth2")])] <- 0
df[c("meth3")][is.na(df[c("meth3")])] <- 0
df[c("repli1")][is.na(df[c("repli1")])] <- mean(df$repli1, na.rm = TRUE)
df[c("repli2")][is.na(df[c("repli2")])] <- mean(df$repli2, na.rm = TRUE)
df[c("GC1k")][is.na(df[c("GC1k")])] <- mean(df$GC1k, na.rm = TRUE)
df[c("recomb")][is.na(df[c("recomb")])] <- mean(df$recomb, na.rm = TRUE)
df[c("atac")][is.na(df[c("atac")])] <- mean(df$atac, na.rm = TRUE)
df[c("CpG_I")][is.na(df[c("CpG_I")])] <- 0
df[c("h3k9me3")][is.na(df[c("h3k9me3")])] <- 0
df[c("h3k36me3")][is.na(df[c("h3k36me3")])] <- 0
##when the downsampling is done, sometimes postions that cant be annotated in certain genomic features occur, to overcome this we just avg them out. 

#transformation of features
df[c("log_recomb")] <- log(df[c("recomb")]+1)

##Loading in model from training
load(model)

#loading in levels
load(levelsval)

new_df <- df[c("context","repli1","GC1k","recomb","meth1","CpG_I","h3k36me3","h3k9me3", "atac")] # make this a input/parameter

# x <- sparse.model.matrix( ~ factor(context, levels = sort(unique(levelsval))) + repli + GC_1k + recomb_decode + meth + CpG_I -1, new_df)
options(na.action='na.pass')

if (log_model == "standard") {
x <- sparse.model.matrix( ~ factor(context, levels = sort(unique(levelsval))) + repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I -1, df) # standard. no interactions
} else if (log_model== "fullinteraction") {
x <- sparse.model.matrix( ~ factor(context, levels = sort(unique(levelsval))) * (repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I) + 
                        (repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I)^2 - 1, df) # interactions between everything
} else if (log_model== "contextinteraction") {
x <- sparse.model.matrix( ~ factor(context, levels = sort(unique(levelsval))) * (repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I) - 1, df) # no interaction between genomic features
}


res_min <- predict(cv.fit, newx = x, s = "lambda.min", type = "response")
res_1se <- predict(cv.fit, newx = x, s = "lambda.1se", type = "response")


#vec1 <- as.vector(unique(df$context))
#vec2 <- as.vector(unique(levelsval))
#unique(vec2[! vec2 %in% vec1])

#predictions_df <-  df
df$prob_min <- res_min
df$prob_1se <- res_1se

write.table(df, file = output, sep='\t', quote=FALSE, row.names = FALSE)