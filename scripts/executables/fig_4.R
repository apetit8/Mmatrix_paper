source("scripts/functions_R/All_functions.R")
library(png)
library(igraph)
library(ggstatsplot)
library(tidyverse)
#####################
sims.dirs <-  list.dirs("simul/fig_4", recursive = FALSE)
#####################

df.fig4 <- df.data(sims.dirs, pattern = "simul/fig_4/", variable="netw", file_size=15000, w_of_6=TRUE, network=TRUE)

plt <- ggbetweenstats(
  data = df.fig4,
  x = pop,
  y = corrM,
  centrality.plotting=FALSE,
  plot.type = "box",
  ggtheme = ggplot2::theme_bw()+theme(plot.margin = unit(c(0.1, 0, 3, 0), "cm")),
  ggplot.component = scale_x_discrete(labels = c('Full Network','Direct regulation\nonly','Distance: 1 node','Distance: 2 nodes','Distance: 3 nodes')),
  pairwise.comparisons=FALSE,
  bf.message=FALSE,
  results.subtitle=FALSE,
  xlab=" ", ylab="Mutational correlation r(M)",
  title = NULL
)
plt


png.netw1 = readPNG('templates/fig_4/networks_fig4.png')

cairo_pdf("figures/fig_4.pdf", width=7, height=5.8)
grid.arrange(
  plt,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
grid.raster(png.netw1, x=0.54, y=0.13, width=0.95)
dev.off()



pfig <- ggplot(data=df.fig4, aes(corrS, corrM))+
  geom_point(aes(col=as.factor(pop)), alpha=0.3, show.legend = FALSE)+
  labs(y="r(M)", x="r(S)", title="Gene number" )+
  coord_fixed(ratio = 1, xlim = c(-1,1), ylim = c(-1,1), expand = TRUE, clip = "on")+
  facet_wrap(pop ~., ncol=6) + theme_bw() #, strip.text = element_blank()
pfig
