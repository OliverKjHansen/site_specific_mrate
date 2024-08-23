#!/usr/bin/env Rscript
library(tidyverse)
library(lmerTest)
library(stringi)
library(glmnet)
#set.seed(14)
"""
This scripts trains a logistic regression model on snp data

As input it takes a file with mutations and path to output

It creates a logistic regression model and a sperate file with the levelorder of the variables.
"""
# inputfiles
args = commandArgs(trailingOnly=TRUE)
#filelists <- list(args[1])
file <- args[1]
output_path <- args[2]
#mutationtype <- args[3]


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Missing arguments (input file).n", call.=FALSE)
}


df <- read_table(file = "data/subsample.txt") %>% mutate(ref=substr(context,(nchar(context)+1)/2,(nchar(context)+1)/2))


##finding the reference and alternative allelse from the data
REF_allel <- unique(df$ref)

nucleotides <- c("A","C","G","T")

find_mutation_type <- function(df, column_name){
  alt_allel <- character()
  for (column in column_name){
      column_data <- df[[column]]
      if (sum(column_data, na.rm = TRUE) > 0) {
      alt_allel <- c(alt_allel, column)
      }
      }
  return(alt_allel)
    }

ALT_allel <- find_mutation_type(df, nucleotides)

# some of the genomic features will have NA values. e.g methylation will not be on a A2N mutation type. I will feed this to the model as a 0. for replication time, C within 1k context, ecombination rate and CpG Island a NA values will just be given the mean of the dataset.
df[c("meth")][is.na(df[c("meth")])] <- 0
df[c("repli")][is.na(df[c("repli")])] <- mean(df$repli, na.rm = TRUE)
df[c("GC_1k")][is.na(df[c("GC_1k")])] <- mean(df$GC_1k, na.rm = TRUE)
df[c("recomb_decode")][is.na(df[c("recomb_decode")])] <- mean(df$recomb_decode, na.rm = TRUE)
df[c("CpG_I")][is.na(df[c("CpG_I")])] <- mean(df$CpG_I, na.rm = TRUE)

#previuos test have shown it is better to log-transform the recombination rate
# should i +1 so i dont get the. -inf after log transforming?
df[c("log_recomb")] <- log(df[c("recomb_decode")]+1)

#before training the model we want to the the context feature, which is closest to the average mutation rate and make it the reference(intercept) for the regression. 
# here is say mutation/opportunities(nonmut+mut)
reference_context <- df %>% 
  group_by(context) %>% 
  summarise(mutations = sum(mut),
         opportunities = n(),
         rate = mutations/opportunities,
         diff = abs(rate-(sum(df[["mut"]])/length(df[["mut"]])))
         ) %>% slice_min(diff, n = 1) %>% as.vector() %>% .$context


## preparing the model

#levelsval <- relevel(factor(df$context),reference_context)
#chaning refernce context
levelsval <- df$context %>% factor() %>% relevel(reference_context)
levelsfile <- paste0(output_path,REF_allel,2,ALT_allel,"_levels.RData")
save(levelsval, file = levelsfile)

y<-df$mut
#x<-data.matrix(df[,c('context','repli','GC_1k','recomb_decode','meth','CpG_I','CpG')])
# x <- model.matrix( ~ context + repli + GC_1k + recomb_decode + meth + CpG_I - 1, df) used to have a -1 to remove intercept 
x <- model.matrix( ~ context + repli + GC_1k + recomb_decode + meth + CpG_I -1, df)
# im not sure what happens to the refernce level when i remove the intercept, but i have to because glmner makes it won intercept read documentation


#grid <- exp(seq(-10.5,-5.5,length=250)) #lambda search space
# C2T keeps on crashing   with my own grid         
#cv.fit <-cv.glmnet(x,y,alpha=1,family='binomial', type.measure = 'deviance', lambda = grid, maxit = 100000) #execute lasso regression

cv.fit <-cv.glmnet(x,y,alpha=1,family='binomial', type.measure = 'deviance', maxit = 1000) # change to 10000 on the cluster

#best_model <- cv.fit

print(cv.fit$lambda.min)
print(cv.fit$lambda.1se)
print(coef(cv.fit, s = "lambda.min"))
print(coef(cv.fit, s = "lambda.1se"))

plot(cv.fit)

filename <- paste0(output_path,REF_allel,2,ALT_allel,"_LassoBestModel.RData")
save(best_model, file = filename) #saving model as an R-object

## these are not necerasary 
#d[c("T")][is.na(d[c("T")])] <- 0
#d[c("A")][is.na(d[c("A")])] <- 0
#d[c("C")][is.na(d[c("C")])] <- 0
#d[c("G")][is.na(d[c("G")])] <- 0
