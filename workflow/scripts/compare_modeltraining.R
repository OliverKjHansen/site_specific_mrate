#!/usr/bin/env Rscript
library("tidyverse")
library("glmnet")

#set.seed(14)
# """
# This scripts trains a logistic regression model on snp data

# As input it takes a file with mutations and path to output

# It creates a logistic regression model and a sperate file with the levelorder of the variables.
# """
# inputfiles
args = commandArgs(trailingOnly=TRUE)
#filelists <- list(args[1])
file <- args[1]
mutationtype <- args[2] 
partition <- args[3]
compare <- args[4]
output <- args[5]


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Missing arguments (input file).n", call.=FALSE)
}

df <- read_table(file = file)

papa <- read_delim(file = partition )
context_dict <- setNames(papa$p_rate, papa$`#pattern`)

# some of the genomic features will have NA values. e.g methylation will not be on a A2N mutation type. I will feed this to the model as a 0. for replication time, C within 1k context, ecombination rate and CpG Island a NA values will just be given the mean of the dataset.
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


#previuos test have shown it is better to log-transform the recombination rate
# should i +1 so i dont get the. -inf after log transforming?
df[c("log_recomb")] <- log(df[c("recomb")]+1)


if (compare == "standard") {
reference_context <- df %>% 
  group_by(context) %>% 
  summarise(mutations = sum(mut),
         opportunities = n(),
         rate = mutations/opportunities,
         diff = abs(rate-(sum(df[["mut"]])/length(df[["mut"]])))
         ) %>% slice_min(diff, n = 1) %>% .[[1,1]]
levelsval <- df$context %>% factor() %>% relevel(reference_context)
df$context <- df$context %>% factor() %>% relevel(reference_context)
y<-df$mut
x <- sparse.model.matrix( ~ context + repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I -1, df)
cv.fit <-cv.glmnet(x,y,alpha=1,family='binomial', type.measure = 'deviance', maxit = 100000)
 # standard. no interactions
} else if (compare == "context") {
reference_context <- df %>% 
  group_by(context) %>% 
  summarise(mutations = sum(mut),
         opportunities = n(),
         rate = mutations/opportunities,
         diff = abs(rate-(sum(df[["mut"]])/length(df[["mut"]])))
         ) %>% slice_min(diff, n = 1) %>% .[[1,1]]
levelsval <- df$context %>% factor() %>% relevel(reference_context)
df$context <- df$context %>% factor() %>% relevel(reference_context)
y<-df$mut
x <- sparse.model.matrix( ~ context -1, df)
cv.fit <-cv.glmnet(x,y,alpha=1,family='binomial', type.measure = 'deviance', maxit = 100000)
} else if (compare == "rate") {
df$context_mapped <- context_dict[df$context]
scalingfactor <- mean(df$mut)/mean(df$context_mapped)
#df$context_mapped <- df$context_mapped * scalingfactor
#df$context_mapped <- df$context_mapped * scalingfactor
df$context_mapped <- pmin(df$context_mapped * scalingfactor, 0.999)
x <- sparse.model.matrix( ~ context_mapped + repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I -1, df)
y<-df$mut
cv.fit <-cv.glmnet(x,y,alpha=1,family='binomial', type.measure = 'deviance', maxit = 100000)
levelsval <- "helloworld"
print(scalingfactor)
save(scalingfactor, file= sub("LassoBestModel", "scalingfactor", output))
} else if (compare == "rateoffset") {
df$context_mapped <- context_dict[df$context]
scalingfactor <- mean(df$mut)/mean(df$context_mapped)
#df$context_mapped <- df$context_mapped * scalingfactor
df$context_mapped <- pmin(df$context_mapped * scalingfactor, 0.999)
x <- sparse.model.matrix( ~ context_mapped + repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I -1, df)
y<-df$mut
cv.fit <-cv.glmnet(x[, -1],y,alpha=1,family='binomial', type.measure = 'deviance', maxit = 100000, offset = x[,"context_mapped"])
levelsval <- "helloworld"
print(scalingfactor)
save(scalingfactor, file= sub("LassoBestModel", "scalingfactor", output))
}

#levelsfile <- sub("models", "levels", output)
levelsfile <- sub("LassoBestModel", "levels", output)
save(levelsval, file = levelsfile)

save(cv.fit, file = output) #saving model as an R-object

