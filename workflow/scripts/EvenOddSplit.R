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
         ) %>% slice_min(diff, n = 1) %>% .[[1,1]]

#changing refernce context
df_even$context <- df_even$context %>% factor() %>% relevel(reference_context[1]) #sometimes there is multiple context with the same rate
df_odd$context <- df_odd$context %>% factor() %>% relevel(reference_context[1])

y_train <- df_even$mut # results for glmnet

if (modeltype == "context") {
x_model <- sparse.model.matrix( ~ context -1, df_even)
y_model <- sparse.model.matrix( ~ context -1, df_odd) #kmer model
} else if (modeltype == "genomic") {
x_model <- sparse.model.matrix( ~ repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac -1, df_even)
y_model <- sparse.model.matrix( ~ repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac -1, df_odd) # genomic model
} else if (modeltype == "linear") {
x_model <- sparse.model.matrix( ~ context + repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac -1, df_even)
y_model <- sparse.model.matrix( ~ context + repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac -1, df_odd) # 
} else if (modeltype == "contextinteraction") {
x_model <- sparse.model.matrix( ~ context * (repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac) -1, df_even)
y_model <- sparse.model.matrix( ~ context * (repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac) -1, df_odd)  # 
} else if (modeltype == "fullinteraction") {
x_model <- sparse.model.matrix( ~ context * (repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac) +
                                      (repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac)^2  -1, df_even)
y_model <- sparse.model.matrix( ~ context * (repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac) +
                                      (repli1 + GC1k + log_recomb + meth1 + CpG_I + h3k9me3 + h3k36me3 + atac)^2  -1, df_odd) # 
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


#save(cv.fit, file = output) #saving model as an R-object

#####PREDICITNG#########

res <- predict(cv.fit , newx = y_model, s = "lambda.min", type = "response")
res_1se <- predict(cv.fit , newx = y_model, s = "lambda.1se", type = "response")
predictions_df <- df_odd
predictions_df$res <- res
predictions_df$res_1se <- res_1se

summary <- predictions_df %>% summarise(obs = n(),
                             k = ncol(x_model)+2, # number of variables and intecept and variance!
                             loglik = sum(log(res)*mut + log(1-res)*(1-mut)),
                             loglik_1se = sum(log(res_1se)*mut + log(1-res_1se)*(1-mut)),
                             AIC = 2*k-2*loglik,
                             AIC_1se = 2*k-2*loglik_1se,
                             BIC = k*log(obs)-2*loglik,
                             BIC_1se = k*log(obs)-2*loglik_1se
                             model =  modeltype,
                             mutationtype = mutationtype)

save(cv.fit, file = gsub("summary", "LassoBestModel", output))
#saving model as an R-object
#save(predictions_df, file = gsub("summary", "predictions", output)) 
#save(summary, file = output)

write.table(summary, file = output, sep='\t', quote=FALSE, row.names = FALSE)




