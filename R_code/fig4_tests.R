# setwd("/home/radillo/Git/GitHub/work/param-infer-clicks/R_code/")
library(ggplot2)
fig4_data <- read.csv(file="../data/fig4_data.csv", header=TRUE, sep=",")

gg <- ggplot(fig4_data,aes(x=trial_nb,y=error)) +
             geom_point(aes(col=model_pair, shape=method),size=6) +
             labs(title="Error as fcn of block size", 
                  y="Relative Error", 
                  x="Block Size (in trials)") +
             scale_colour_brewer(palette = "Set1") # change color palette
#gg + facet_wrap(method~model_pair, nrow = 2)
plot(gg + theme(text = element_text(size=30)) + 
  facet_grid(method~model_pair) + 
    scale_x_discrete(expand=c(0.1, 0.1),
                     labels=c("","100","300","","500"),
                     limits=unique(fig4_data$trial_nb)))
  
#gg + theme_bw()  # change whole theme
#gg)e