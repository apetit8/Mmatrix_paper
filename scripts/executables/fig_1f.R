source("scripts/functions_R/All_functions.R")
#########################################
angle = c(0, 0.393, 0.785, 1.178)
#####################

M.factor <- 1/2
S.factor <- 1/1500


cairo_pdf("figures/fig1_f.pdf", width=3.5, height=3.6)
laymat <- matrix(seq(1,16, by=1), nrow = 4, ncol = 4)
par(mar=c(0.1,0.1,0.1,0.1), oma=c(0, 0, 0, 0))
layout(mat = laymat) # Widths of the two columns

for (i in angle) {
  oneplot.allellipse(sprintf("simul/fig_1abc/m/simuangle%s",i), S.factor=S.factor, M.factor=M.factor, Sell=TRUE, Mell=FALSE,
                   xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), all.reps=FALSE, xlab="", ylab="",
                   yaxt = "n", xaxt = "n", legend=FALSE, asp=1, axes = FALSE,mgp = c(0, 0, 0))
  if(i==0){text("Selection", x=0, y=0.15)}
  oneplot.allellipse(sprintf("simul/fig_1abc/m/simuangle%s",i), S.factor=S.factor, M.factor=M.factor, Sell=FALSE, Mell=TRUE,
                   xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), all.reps=FALSE, xlab="", ylab="",
                   yaxt = "n", xaxt = "n", legend=FALSE, asp=1, axes = FALSE,mgp = c(0, 0, 0))
  if(i==0){text("M multilinear", x=0, y=0.15)}
  oneplot.allellipse(sprintf("simul/fig_1abc/fkl/simuangle%s",i), S.factor=S.factor, M.factor=M.factor, Sell=FALSE, Mell=TRUE,
                     xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), all.reps=FALSE, xlab="", ylab="",
                     yaxt = "n", xaxt = "n", legend=FALSE, asp=1, axes = FALSE,mgp = c(0, 0, 0))
  if(i==0){text("M FKL", x=0, y=0.15)}
  oneplot.allellipse(sprintf("simul/fig_1abc/w/simuangle%s",i), S.factor=S.factor, M.factor=M.factor, Sell=FALSE, Mell=TRUE,
                   xlim=c(-0.2,0.2), ylim=c(-0.2,0.2), all.reps=FALSE, xlab="", ylab="",
                   yaxt = "n", xaxt = "n", legend=FALSE, asp=1, axes = FALSE,mgp = c(0, 0, 0))
  if(i==0){text("M GRN", x=0, y=0.15)}
}
dev.off()




