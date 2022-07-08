source("../functions_R/All_functions.R")
library(png)
library(igraph)
library(ggstatsplot)
library(tidyverse)

#####################
sims.dirs1 <-  c("../../simul/fig_4/indirect-1","../../simul/fig_4/indirect-2","../../simul/fig_4/indirect-3")
sims.dirs2 <-  list.dirs("../../simul/fig_supp2", recursive = FALSE)
#####################
d1 <- df.data(sims.dirs1, pattern = "../../simul/fig_4/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE)
d2 <- df.data(sims.dirs2, pattern = "../../simul/fig_supp2/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE)
df.figs2 <- rbind(d1,d2)

plt <- ggbetweenstats(
  data = df.figs2,
  x = pop,
  y = corr,
  centrality.plotting=FALSE,
  plot.type = "box",
  ggtheme = ggplot2::theme_bw()+theme(plot.margin = unit(c(0.1, 0, 3, 0), "cm")),
  # ggplot.component = scale_x_discrete(labels = c('Full Network','Direct only','1 intermediate','2 intermediate','3 intermediate')),
  pairwise.comparisons=FALSE,
  bf.message=FALSE,
  results.subtitle=FALSE,
  xlab=" ", ylab="Mutational correlation (rM)",
  title = NULL
)
plt


png.netw1 = readPNG('../../templates/fig_4/networks_figs2.png')

cairo_pdf("../../figures/fig_supp2.pdf", width=7, height=5.8)
grid.arrange(
  plt,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
grid.raster(png.netw1, x=0.54, y=0.13, width=0.95)
dev.off()

