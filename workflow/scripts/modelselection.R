library("tidyverse")

args = commandArgs(trailingOnly=TRUE)


#looping through inout and make one big file

modelselection <- tibble()
for (file in args) { #minus one because the last is path to output
  load(file)
  modelselection <- rbind(modelselection,file)
}

pbic <- modelselection %>%  ggplot(mapping = aes(x = model, y = AIC)) + 
  theme_pubr() + 
  facet_wrap(.~mutationtype, scales = "free") +
  geom_point()


pbic <- modelselection %>%  ggplot(mapping = aes(x = model, y = BIC)) + 
  theme_pubr() + 
  facet_wrap(.~mutationtype, scales = "free") +
  geom_point()