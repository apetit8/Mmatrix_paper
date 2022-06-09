source("../functions_R/All_functions.R")
library(png)
#########################################
sims.dirs <- list.dirs("../../simul/fig_4", recursive = FALSE)

angle = c(0.8, 0, -0.4)
of        <- "fig1"
modulo <- pi
lim <- c(-0.2,0.2)
#####################

M.factor <- 1
S.factor <- 0.0008

# cairo_pdf("../../figures/fig_4_ellipses.pdf", width=10, height=8)
cairo_pdf("../../figures/fig_4_ellipses.pdf", width=7, height=6.1)
laymat <- matrix(seq(1,15, by=1), nrow = 3, ncol = 5)
par(mar=c(0.1,0.1,0.1,0.1), oma=c(10, 0, 0, 0))
layout(mat = laymat) # Widths of the two columns
for (sims.dir in sims.dirs) {
  oneplot.netw(sims.dir, G.factor=G.factor, M.factor=M.factor, S.factor=S.factor,
                     # main= paste0("\u03BE", "\u03B1", sprintf(" : m = %s, w = %s", mse.m, mse.w)),
                     xlim=lim, ylim=lim, all.reps=TRUE, all.gen=FALSE,
                     yaxt = "n", xaxt = "n", mgp = c(1, 1, 0), asp=1, angle = angle)
}

png.netw1 = readPNG('../../templates/fig_4/netw_topo.png')
grid.raster(png.netw1, x=0.5, y=0.085, width=0.95)

dev.off()
print("Ellipses done !")


