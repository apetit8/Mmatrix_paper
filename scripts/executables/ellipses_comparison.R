source("../functions_R/figure_tools.R")

#####################
sims.dirs <- c("../../simul/fig_2/c0_0adda",
               "../../simul/fig_2/round_s"
               )
#angle = c(-0.2, 1.2)
angle = c(0.8) #c(-1.4, -1.3, -1.2 , -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4)
of        <- "fig1"
where     <- "figures"
#####################


M.factor <- 4
G.factor <- 0.25
#S.factor <- 0.022 #Multilinear
#S.factor <- 0.01 #Multilinear scaled 0_1
#S.factor <- 0.0008 #Wagner
S.factor <- 0.0044

# pdfname   <- print(sprintf("../../figures/%s_ellipse.pdf", of))
# #Draw ellipses
# pdf(pdfname, width=13, height=4.4)
# layout(t(1:3))

j <- 1
for (i in angle) {

  sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle_%s.par", i), full.names=TRUE)
  
  # oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor,asp=1, main=sprintf("MSE m = %s MSE w = %s", mse.m, mse.w),xlim=c(-0.5,0.5), ylim=c(-0.5,0.5), all.reps=TRUE, Gell=FALSE) #xlim=c(-3,3), ylim=c(-3,3) all.gen=TRUE all.reps=TRUE
  oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor,asp=1, main=(i),xlim=c(-0.5,0.5), ylim=c(-0.5,0.5), all.reps=FALSE, Gell=FALSE)
}
# dev.off()
# print("Ellipses done !")
