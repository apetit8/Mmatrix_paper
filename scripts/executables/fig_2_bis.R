source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs.fig2 <-list.dirs(
  "../../simul/fig_2_bis", recursive = FALSE
)
modulo <- pi
#####################

#Data
df.fig2 <- df.data(sims.dirs.fig2, pattern = "../../simul/fig_2/d", variable="netw")
g_sdfig2 <- table.index(sims.dirs.fig2, pattern = "../../simul/fig_2/", ref = df.fig2, bymean = FALSE, asgrob=FALSE)

#Plot
pfig2 <- ggplot(data=g_sdfig2, aes(Population, Xi_alpha, fill=Population))+
  geom_col()
pfig2
###############PDF
cairo_pdf("../../figures/fig_2.pdf", width=10, height=8)

grid.arrange(
  pfig2,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
grid.raster(png.netw1, x=.1, y=0.75, width=0.14)
grid.raster(png.netw2, x=.93, y=0.75, width=0.14)
grid.raster(png.netw3, x=.1, y=0.3, width=0.14)
grid.raster(png.netw4, x=.93, y=0.3, width=0.14)

dev.off()


cairo_pdf("../../figures/fig_2_small.pdf", width=7.5, height=6)

grid.arrange(
  pfig2,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
grid.raster(png.netw1, x=.1, y=0.75, width=0.14)
grid.raster(png.netw2, x=.92, y=0.75, width=0.14)
grid.raster(png.netw3, x=.1, y=0.3, width=0.14)
grid.raster(png.netw4, x=.92, y=0.3, width=0.14)

dev.off()
