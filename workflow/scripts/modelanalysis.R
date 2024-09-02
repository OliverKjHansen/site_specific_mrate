"""
This script should take logistic regression models and make a preformance analysis 
"""

library("tidyverse")
library("glmnet")
library("ggpubr")
library("broom")
library("plotmo")
library("ggrepel")


#make this de
if snp =
mutationtypes <- c("C2A","C2T","C2G","A2T","A2C","A2G")

if indel =
mutationtypes <- c("insertion","deletion")
#g_features <- c("repli","GC_1k","recomb_decode","CpG_I","meth")

tidied_cv <- tibble()
glanced_cv <- tibble()
for (type in mutationtypes) {
  df <- load(paste0("data/models/",type,"_LassoBestModel.RData"))
  glanced <- glance(cv.fit) %>% mutate(mutationtype = type)
  tibble <- tidy(cv.fit) %>% mutate(mutationtype = type)
  glanced_cv <- rbind(glanced_cv,glanced)
  tidied_cv <- rbind(tidied_cv,tibble)
}


tidied_cv %>% ggplot(aes(lambda, estimate)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), colour="grey") +
  geom_point(colour="red") +
  facet_wrap( .~mutationtype, nrow = 2 , ncol = 3,scales = "free") +
  geom_vline(data = glanced_cv, aes(xintercept = lambda.min), linetype = "dotted", colour = "black") +
  geom_vline(data = glanced_cv, aes(xintercept = lambda.1se), linetype = "dotted", colour = "black") +
  xlab("Lambda (λ)") +
  ylab("Crossvalidated Binomial Deviance") +
  scale_x_log10() +
  theme_classic()


tidied_cv %>% ggplot(aes(lambda, nzero)) +
  geom_line() +
  facet_wrap( .~mutationtype, nrow = 2 , ncol = 3,scales = "free") +
  geom_vline(data = glanced_cv, aes(xintercept = lambda.min), linetype = "dotted", colour = "black") +
  geom_vline(data = glanced_cv, aes(xintercept = lambda.1se), linetype = "dotted", colour = "black") +
  xlab("Lambda (λ)") +
  ylab("#Not Zero Coefficients") +
  scale_x_log10() +
  theme_classic()

tidied_coef <- tibble()
#glanced_coef <- tibble()
for (type in mutationtypes) {
  df <- load(paste0("data/models/",type,"_LassoBestModel.RData")) # change with input
  #gl_coef <- glance(cv.fit$glmnet.fit) %>% mutate(mutationtype = type)
  tid_coef <- tidy(cv.fit$glmnet.fit,return_zeros = TRUE) %>% mutate(mutationtype = type)
  #glanced_coef <- rbind(glanced_coef,glanced_coef)
  tidied_coef <- rbind(tidied_coef,tid_coef)
}


tidied_coef %>%
  filter(term != "(Intercept)") %>%
  mutate(term = gsub("context","", term)) %>% 
  group_by(mutationtype) %>% 
  mutate(minimum = if_else(lambda == min(lambda), TRUE, FALSE)) %>% 
  ungroup() %>% 
  group_by(mutationtype,lambda) %>% 
  #mutate(outlier = if_else((minimum == TRUE & (estimate < quantile(estimate,0.25)-1.2*IQR(estimate) | estimate > quantile(estimate,0.75)+1.2*IQR(estimate))), term, NA)) %>% 
  mutate(outlier = if_else(minimum == TRUE & (rank(-estimate) %in% 1:3 | rank(estimate)%in% 1:3), term, NA)) %>%
  ggplot(aes(lambda, estimate, group = term)) +
  geom_text_repel(aes(label = outlier), colour = "red", na.rm = TRUE) + 
  geom_line() +
  facet_wrap( .~mutationtype, nrow = 2 , ncol = 3,scales = "free") +
  scale_x_log10() +
  theme_classic() +
  xlab("Lambda (λ)") +
  ylab("Coefficient")


#not finished
# models_df %>% 
#   mutate(coefficient = ifelse(estimate > 0, "Positive","Negative")) %>% 
#   filter(criteria == "lambda.1se",
#          term %in% g_features) %>%
#   ggplot(mapping = aes(x=term, y=mutationtype, fill = coefficient)) + geom_tile() +  
#   scale_fill_manual(values = c("Positive" = "red", "Negative" = "darkblue")) +
#   theme_pubr()

