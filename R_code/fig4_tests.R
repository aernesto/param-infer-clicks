# setwd("/home/radillo/Git/GitHub/work/param-infer-clicks/R_code/")
library(ggplot2)
fig4_data_raw <- read.csv(file="../data/fig4_data.csv", header=TRUE, sep=",")
fig4_data <- subset(fig4_data_raw, method == "stoch")
fig4_data$model_pair = factor(
  fig4_data$model_pair, levels = c("linlin","linnonlin","nonlinlin","nonlinnonlin"),
  labels = c("L-L", "L-NL", "NL-L","NL-NL")
)
gg <- ggplot(fig4_data,aes(x=trial_nb,y=error)) +
             geom_line(aes(col=model_pair),size=2) +
             geom_point(aes(col=model_pair), size=4) +
             labs(y="Relative Error", 
                  x="Block Size (in trials)") +
             scale_colour_brewer(palette = "Set1") # change color palette

#gg
#gg + facet_wrap(method~model_pair, nrow = 2)

# plot(gg + 
#   facet_grid(method~model_pair) + 
#     scale_x_discrete(expand=c(0.16, 0.16),
#                      labels=c("","100","300","","500"),
#                      limits=unique(fig4_data$trial_nb)) +
#     theme(legend.position="top"))
plot(gg+theme(text = element_text(size=30)))
      #gg + theme_bw()  # change whole theme
#gg)

