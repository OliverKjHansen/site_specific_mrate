
library(tidyverse)
#library(stringi)
library(glmnet)

# inputfile
args = commandArgs(trailingOnly=TRUE)

file <- args[1]
model <- args[2]
levelsval <- args[3]
output <- args[4]

# loading in prediction dataset
df <- read.table(file, header = TRUE)
df[c("meth")][is.na(df[c("meth")])] <- 0
df[c("repli")][is.na(df[c("repli")])] <- mean(df$repli, na.rm = TRUE)
df[c("GC_1k")][is.na(df[c("GC_1k")])] <- mean(df$GC_1k, na.rm = TRUE)
df[c("recomb_decode")][is.na(df[c("recomb_decode")])] <- mean(df$recomb_decode, na.rm = TRUE)
df[c("CpG_I")][is.na(df[c("CpG_I")])] <- mean(df$CpG_I, na.rm = TRUE)
##when the downsampling is done, sometimes postions that cant be annotated in certain genomic features occur, to overcome this we just avg them out. 

#transformation of features
df[c("recomb_decode")] <- log(df[c("recomb_decode")]+1)

##Loading in model from training
load(model)

#loading in levels
load(levelsval)

new_df <- df[c("context","repli","GC_1k","recomb_decode","meth","CpG_I")] # make this a input/parameter

x <- model.matrix( ~ factor(context, levels = sort(unique(levelsval))) + repli + GC_1k + recomb_decode + meth + CpG_I -1, new_df)
options(na.action='na.pass')

# if (log_model == "nobeta") {
# x <- model.matrix( ~ context + repli + GC_1k + recomb_decode + meth + CpG_I -1, df) # no intercepts
# } else if (log_model== "intercept") {
# x <- model.matrix( ~ context + repli + GC_1k + recomb_decode + meth + CpG_I, df) # with intercept
# }


res_min <- predict(cv.fit, newx = x, s = "lambda.min", type = "response")
res_1se <- predict(cv.fit, newx = x, s = "lambda.1se", type = "response")


#vec1 <- as.vector(unique(df$context))
#vec2 <- as.vector(unique(levelsval))
#unique(vec2[! vec2 %in% vec1])

predictions_df <-  df
predictions_df$prob_min <- res_min
predictions_df$prob_1se <- res_1se

write.table(predictions_df, file = output, sep='\t', quote=FALSE)