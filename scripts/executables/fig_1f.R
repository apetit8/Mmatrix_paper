source("../functions_R/All_functions.R")
#########################################
sims.dirs <- list.dirs("../../simul/fig_1", recursive = FALSE)
angle = c(0, 0.4, 0.8, 1.2)
of        <- "fig1"
modulo <- pi
#####################

M.factor <- 1500
G.factor <- 0.25
S.factor <- 1

laymat <- matrix(seq(1,12, by=1), nrow = 3, ncol = 4)
par(mar=c(0.1,0.1,0.1,0.1), oma=c(10, 0, 0, 0))
layout(mat = laymat) # Widths of the two columns

for (i in angle) {
oneplot.allellipse(sprintf("../../simul/fig_1/m/simuangle%s.parangle%s.par", i,i), G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, Sell=TRUE, Mell=FALSE,
                   xlim=c(-7,7), ylim=c(-7,7), all.reps=FALSE, xlab="", ylab="",
                   yaxt = "n", xaxt = "n", legend=legend, asp=1, axes = FALSE,mgp = c(0, 0, 0))
oneplot.allellipse(sprintf("../../simul/fig_1/m/simuangle%s.parangle%s.par", i,i), G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, Sell=FALSE, Mell=TRUE,
                   xlim=c(-7,7), ylim=c(-7,7), all.reps=FALSE, xlab="", ylab="",
                   yaxt = "n", xaxt = "n", legend=legend, asp=1, axes = FALSE,mgp = c(0, 0, 0))
oneplot.allellipse(sprintf("../../simul/fig_1/w/simuangle%s.parangle%s.par", i,i), G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, Sell=FALSE, Mell=TRUE,
                   xlim=c(-7,7), ylim=c(-7,7), all.reps=FALSE, xlab="", ylab="",
                   yaxt = "n", xaxt = "n", legend=legend, asp=1, axes = FALSE,mgp = c(0, 0, 0))
}
