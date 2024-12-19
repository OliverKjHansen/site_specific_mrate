library(tidyverse)
#library(stringi)
library(glmnet)

# inputfile
args = commandArgs(trailingOnly=TRUE)

file <- args[1]
partition <- args[2]
model <- args[3]
compare <-args[4]
levelsval <- args[5]
output <- args[6]

papa <- read_delim(file = partition)
context_dict <- setNames(papa$p_rate, papa$`#pattern`)

#load(sub("LassoBestModel", "scalingfactor", model))

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
#load(model)

#loading in levels
#load(levelsval)

#new_df <- df[c("context","repli1","GC1k","recomb","meth1","CpG_I","h3k36me3","h3k9me3", "atac")] # make this a input/parameter

# x <- sparse.model.matrix( ~ factor(context, levels = sort(unique(levelsval))) + repli + GC_1k + recomb_decode + meth + CpG_I -1, new_df)
#options(na.action='na.pass')

if (compare == "standard") {
load(model)
load(levelsval)
#y<-df$mut
x <- sparse.model.matrix( ~ factor(context, levels = sort(unique(levelsval))) + repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I -1, df) # standard. no interactions
res_model <- predict(cv.fit, newx = x, s = "lambda.min", type = "response")
df$context_mapped <- context_dict[df$context] # this column is actually the predictions from kmerpapa
} else if (compare == "context") {
load(model)
load(levelsval)
#y<-df$mut
x <- sparse.model.matrix( ~ factor(context, levels = sort(unique(levelsval))) -1, df) # standard. no interactions
res_model <- predict(cv.fit, newx = x, s = "lambda.min", type = "response")
df$context_mapped <- context_dict[df$context] # this column is actually the predictions from kmerpapa
} else if (compare == "rate") {
load(sub("LassoBestModel", "scalingfactor", model))
load(model)
df$context_mapped <- context_dict[df$context]
#df$context_mapped <- df$context_mapped * scalingfactor
df$context_mapped <- pmin(df$context_mapped * scalingfactor, 0.999)
x <- sparse.model.matrix( ~ context_mapped + repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I -1, df)
#y<-df$mut
res_model <- predict(cv.fit, newx = x, s = "lambda.min", type = "response")
} else if (compare == "rateoffset") {
load(sub("LassoBestModel", "scalingfactor", model))
load(model)
df$context_mapped <- context_dict[df$context]
#df$context_mapped <- df$context_mapped * scalingfactor
df$context_mapped <- pmin(df$context_mapped * scalingfactor, 0.999)
x <- sparse.model.matrix( ~ context_mapped + repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I -1, df)
#y<-df$mut
#cv.fit <-cv.glmnet(x[, -1],y,alpha=1,family='binomial', type.measure = 'deviance', maxit = 100000, offset = x[,"context_mapped"])
res_model<- predict(cv.fit, newx = x[,2:ncol(x)], newoffset = df$context_mapped %>% as.vector(), s = "lambda.min", type = "response")
}


#res_papa <- predict(cv.fit, newx = x, s = "lambda.min", type = "response")
#res_model <- predict(cv.fit, newx = x, s = "lambda.1se", type = "response")


#vec1 <- as.vector(unique(df$context))
#vec2 <- as.vector(unique(levelsval))
#unique(vec2[! vec2 %in% vec1])

#predictions_df <-  df
df$prob_papa <- df$context_mapped # this is the kmaerpapa predictions
df$prob_model <- res_model

final <- df %>% select(c("chrom","pos", "mut", "prob_papa" ,"prob_model"))

write.table(final, file = output, sep='\t', quote=FALSE, row.names = FALSE)