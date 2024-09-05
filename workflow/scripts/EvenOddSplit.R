#!/usr/bin/env Rscript

#!/usr/bin/env Rscript
library("tidyverse")
library("glmnet")
#set.seed(14)
# inputfile
args = commandArgs(trailingOnly=TRUE)
file <- args[1]
modeltype <- args[2]
mutationtype <-args[3]
output <- args[4]


df <- read_table(file)

# ADDING CERTAIN EDITS TO DATAFRAME
df[c("meth")][is.na(df[c("meth")])] <- 0
df[c("repli")][is.na(df[c("repli")])] <- mean(df$repli, na.rm = TRUE)
df[c("GC_1k")][is.na(df[c("GC_1k")])] <- mean(df$GC_1k, na.rm = TRUE)
df[c("recomb_decode")][is.na(df[c("recomb_decode")])] <- mean(df$recomb_decode, na.rm = TRUE)
df[c("CpG_I")][is.na(df[c("CpG_I")])] <- mean(df$CpG_I, na.rm = TRUE)

#transformation of features
df[c("recomb_decode")] <- log(df[c("recomb_decode")]+1)

# Create a logical vector to check if the numeric part is even
is_even <- as.integer(gsub("chr", "", df$chrom)) %% 2 == 0

# Split the data frame into two based on the condition
df_even <- df[is_even, ] # dataframe with only even chromosomes. We train on this 

df_odd <- df[!is_even, ] # dataframe with only odd chromosomes. We predict on this


reference_context <- df_even %>% 
  group_by(context) %>% 
  summarise(mutations = sum(mut),
         opportunities = n(),
         rate = mutations/opportunities,
         diff = abs(rate-(sum(df[["mut"]])/length(df[["mut"]])))
         ) %>% slice_min(diff, n = 1) %>% as.vector() %>% .$context

#changing refernce context
df_even$context <- df_even$context %>% factor() %>% relevel(reference_context[1]) #sometimes there is multiple context with the same rate
df_odd$context <- df_odd$context %>% factor() %>% relevel(reference_context[1])

y_train <- df_even$mut # results for glmnet

if (modeltype == "kmer") {
x_model <- model.matrix( ~ context -1, df_even)
y_model <- model.matrix( ~ context -1, df_odd) #kmer model
} else if (modeltype == "genomic")  {
x_model <- model.matrix( ~ repli + GC_1k + recomb_decode + meth + CpG_I -1, df_even)
y_model <- model.matrix( ~ repli + GC_1k + recomb_decode + meth + CpG_I -1, df_odd) # genomic model
} else if (modeltype == "both") {
x_model <- model.matrix( ~ context + repli + GC_1k + recomb_decode + meth + CpG_I -1, df_even) 
y_model <- model.matrix( ~ context + repli + GC_1k + recomb_decode + meth + CpG_I -1, df_odd) # complex model
}

# if the odd and even chromo have a different set of context we force the missing in but with all zeros. else the model breaks
not_in_y <- setdiff(colnames(x_model), colnames(y_model))
not_in_x <- setdiff(colnames(y_model), colnames(x_model))

for (col in not_in_y) {
  y_model <- cbind(y_model, col = 0)
}

for (col in not_in_x) {
  x_model <- cbind(x_model, col = 0)
}

cv.fit <- cv.glmnet(x_model, y_train, alpha=1,family='binomial', type.measure = 'deviance', maxit = 100000)


save(cv.fit, file = output) #saving model as an R-object

#####PREDICITNG#########

res <- predict(cv.fit , newx = y_model, s = "lambda.min", type = "response")
predictions_df <- df_odd
predictions_df$res <- res

summary <- predictions_df %>% summarise(obs = n(),
                             k = ncol(x_model)+2, # number of variables and intecept and variance!
                             loglik = sum(log(res)*mut + log(1-res)*(1-mut)),
                             AIC = 2*k-2*loglik,
                             BIC = k*log(obs)-2*loglik,
                             model =  modeltype,
                             mutationtype = mutationtype)

save(cv.fit, file = gsub("summary", "LassoBestModel", output))
 #saving model as an R-object
save(predictions_df, file = gsub("summary", "predictions", output)) 

save(summary, file = output)




