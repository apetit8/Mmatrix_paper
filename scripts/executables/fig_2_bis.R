source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs.fig2 <-list.dirs(
  "../../simul/fig_2_bis", recursive = FALSE
)

modulo <- pi
#####################

#Data
df.fig2 <- df.data(sims.dirs.fig2, pattern = "../../simul/fig_2_bis/d", variable="netw")
g_sdfig2 <- table.index(sims.dirs.fig2, pattern = "../../simul/fig_2_bis/d", ref = df.fig2, bymean = FALSE, asgrob=FALSE)
g_sdfig2[,5] <- c( mean(subset(df.fig2, netw == "0")$ecc_M),
                   mean(subset(df.fig2, netw == "1")$ecc_M),
                mean(subset(df.fig2, netw == "2")$ecc_M),
                mean(subset(df.fig2, netw == "3")$ecc_M),
                mean(subset(df.fig2, netw == "4")$ecc_M),
                mean(subset(df.fig2, netw == "5")$ecc_M))
names(g_sdfig2) <- c("Population","sq_rho","standard deviation","Xi_alpha", "Eccentricity")

#Plot
pfig2 <- ggplot(data=g_sdfig2, aes(Population, Xi_alpha))+
  geom_col()+
  labs(y="\u03BE\u03B1", x="Connectivity distance")+
  theme(plot.margin = unit(c(0.3,8, 0.3, 0.3), "cm"))
pfig2

png.dist.netw = readPNG('../../templates/fig_2_bis/distance_network2.png')

###############PDF
cairo_pdf("../../figures/fig_2_bis.pdf", width=6, height=3)
grid.arrange(
  pfig2,
  ncol = 1,
  nrow = 1,
  clip = FALSE
)

grid.raster(png.dist.netw, x=0.75, y=0.5, width=0.5)

dev.off()


cairo_pdf("../../figures/fig_2_ter.pdf", width=9, height=3)
grid.arrange(
  pfig2b, pfig2,
  ncol = 2,
  nrow = 1,
  widths = c(0.5, 1),
  clip = FALSE
)
grid.raster(png.dist.netw, x=0.84, y=0.5, width=0.33)
dev.off()

