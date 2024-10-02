library("tidyverse")
library("plotly")
library("svglite")
library("ggpubr")
# """""""
# takes a parameter grid file and findes the best set of parameters while also plot the parameter search
# """""""


args = commandArgs(trailingOnly=TRUE)
mutationtype <- args[1]
inputfile <- args[2]
output_plot <-args[3]
output_parameters <- args[4]

df <- read.table(file = inputfile, header = TRUE) %>% 
  mutate(minimum = if_else(LL_test == min(LL_test), TRUE, FALSE),#have to check if it is loglikelihood  or negative logloikelihood
  P = as.factor(P),
  alpha = as.factor(alpha))

# df %>% 
#   ggplot(mapping = aes(x = alpha, y = LL_test, colour = minimum)) +
#   geom_point() + 
#   theme_pubr() + 
#   ggtitle(mutationtype) +
#   theme(plot.title = element_text(hjust = 0.5))
  
line_plot <- df %>% 
  ggplot(mapping = aes(x = P, y = LL_test, colour = alpha)) +
  geom_line(aes(group = alpha)) + geom_point() + 
  theme_pubr() + 
  ggtitle(mutationtype) +
  theme(plot.title = element_text(hjust = 0.5))

grid_plot <- df %>% 
  ggplot(mapping = aes(x = alpha, y = P, colour = minimum)) + #maybe change?
  geom_point() + 
  theme_pubr(legend = "right") + 
  ggtitle(mutationtype) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Alpha (pseudo count)") +
  ylab("Penalty (complexity)") +
  guides(colour=guide_legend(title.position = "top", ncol=1))
  
# third_dimension <- plot_ly(df, x = ~alpha, y = ~P, z = ~LL_test) %>% 
#   add_markers(color = ~LL_test)

best_set <- df %>% filter(minimum == TRUE) %>% select(c("alpha","P"))

#best_alpha <- best_set[[1,"alpha"]]
#best_pseudo <- best_set[[1,"P"]]

write_delim(best_set, file = output_parameters, col_names = FALSE)


ggsave(sub(".pdf", ".svg", output_plot), plot = line_plot) # as svg
ggsave(sub(".pdf", ".jpg", output_plot), plot = line_plot) #as jpeg
ggsave(output_plot, plot = line_plot, width = 10, height = 10) # as pdf
ggsave(output_plot, plot = grid_plot, width = 10, height = 10) # as pdf


#pdf(file = "FileName.pdf", width = 8, height = 11) # defaults to 7 x 7 inches
#plot(x, y)
#dev.off()



