library(ggplot2)
fig4_data <- read.csv(file="../data/fig4_data.csv", header=TRUE, sep=",")

gg <- ggplot(fig4_data,aes(x=trial_nb,y=error)) +
             geom_point(aes(col=model_pair, shape=method),size=4) +
             labs(title="Error as fcn of block size", 
                  subtitle="Iteration 1", y="Relative Error", 
                  x="Block Size (in trials)") +
             scale_colour_brewer(palette = "Set1") # change color palette
#gg + facet_wrap(method~model_pair, nrow = 2)
plot(gg + theme(text = element_text(size=20)) + 
  facet_grid(method~model_pair) + 
    scale_x_discrete(expand=c(0.1, 0.1), 
                     limits=unique(fig4_data$trial_nb)))
  
#gg + theme_bw()  # change whole theme
#gg)e