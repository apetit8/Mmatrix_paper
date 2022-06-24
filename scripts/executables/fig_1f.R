source("../functions_R/All_functions.R")
#########################################
angle = c(0, 0.393, 0.785, 1.178)
#####################

M.factor <- 1
S.factor <- 1/1500

laymat <- matrix(seq(1,12, by=1), nrow = 3, ncol = 4)
par(mar=c(0.1,0.1,0.1,0.1), oma=c(10, 0, 0, 0))
layout(mat = laymat) # Widths of the two columns

cairo_pdf("../../figures/fig1_f.pdf", width=3.5, height=2.6)
laymat <- matrix(seq(1,12, by=1), nrow = 3, ncol = 4)
par(mar=c(0.1,0.1,0.1,0.1), oma=c(0, 0, 0, 0))
layout(mat = laymat) # Widths of the two columns

for (i in angle) {
  oneplot.allellipse(sprintf("../../simul/fig_1/m/simuangle%s.parangle%s.par", i,i), S.factor=S.factor, M.factor=M.factor, Sell=TRUE, Mell=FALSE,
                   xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), all.reps=FALSE, xlab="", ylab="",
                   yaxt = "n", xaxt = "n", legend=legend, asp=1, axes = FALSE,mgp = c(0, 0, 0))
  oneplot.allellipse(sprintf("../../simul/fig_1/m/simuangle%s.parangle%s.par", i,i), S.factor=S.factor, M.factor=M.factor, Sell=FALSE, Mell=TRUE,
                   xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), all.reps=FALSE, xlab="", ylab="",
                   yaxt = "n", xaxt = "n", legend=legend, asp=1, axes = FALSE,mgp = c(0, 0, 0))
  oneplot.allellipse(sprintf("../../simul/fig_1/w/simuangle%s.parangle%s.par", i,i), S.factor=S.factor, M.factor=M.factor, Sell=FALSE, Mell=TRUE,
                   xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), all.reps=FALSE, xlab="", ylab="",
                   yaxt = "n", xaxt = "n", legend=legend, asp=1, axes = FALSE,mgp = c(0, 0, 0))
}
dev.off()




