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
log_model <- args[3] # do we want to make a intercept or not
output <- args[4]


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Missing arguments (input file).n", call.=FALSE)
}

df <- read_table(file = file)

#df <- read_table(file = file) %>% mutate(ref=substr(context,(nchar(context)+1)/2,(nchar(context)+1)/2))

##finding the reference and alternative allelse from the data
#REF_allel <- unique(df$ref)

#nucleotides <- c("A","C","G","T")

# find_mutation_type <- function(df, column_name){
#   alt_allel <- character()
#   for (column in column_name){
#       column_data <- df[[column]]
#       if (sum(column_data, na.rm = TRUE) > 0) {
#       alt_allel <- c(alt_allel, column)
#       }
#       }
#   return(alt_allel)
#     }

#ALT_allel <- find_mutation_type(df, nucleotides)

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

#before training the model we want to the the context feature, which is closest to the average mutation rate and make it the reference(intercept) for the regression. 
# here is say mutation/opportunities(nonmut+mut)
reference_context <- df %>% 
  group_by(context) %>% 
  summarise(mutations = sum(mut),
         opportunities = n(),
         rate = mutations/opportunities,
         diff = abs(rate-(sum(df[["mut"]])/length(df[["mut"]])))
         ) %>% slice_min(diff, n = 1) %>% .[[1,1]]


## preparing the model

#levelsval <- relevel(factor(df$context),reference_context)
#chaning reference context
levelsval <- df$context %>% factor() %>% relevel(reference_context)
df$context <- df$context %>% factor() %>% relevel(reference_context)

levelsfile <- sub("models", "levels", output)
levelsfile <- sub("LassoBestModel", "levels", levelsfile) # making output name

# print(log_model,mutationtype)
# print(levelsfile)

save(levelsval, file = levelsfile)


## lets try to make two models. one where number of dummy variables are k and one where it is k-1. basically do we choose to fit the intercept or not?
y<-df$mut

if (log_model == "standard") {
x <- sparse.model.matrix( ~ context + repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I -1, df) # standard. no interactions
} else if (log_model== "fullinteraction") {
x <- sparse.model.matrix( ~ context * (repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I) + 
                            (repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I)^2 - 1, df) # interactions between everything
} else if (log_model== "contextinteraction") {
x <- sparse.model.matrix( ~ context * (repli1 + GC1k + log_recomb + meth1 + atac + h3k9me3 + h3k36me3 + CpG_I) - 1, df) # no interaction between genomic features
}

print(colnames(x))
#f <- as.formula( ~ .*.) and x <- model.matrix(f, TrainData)[, -1]

cv.fit <-cv.glmnet(x,y,alpha=1,family='binomial', type.measure = 'deviance', maxit = 100000)

save(cv.fit, file = output) #saving model as an R-object

