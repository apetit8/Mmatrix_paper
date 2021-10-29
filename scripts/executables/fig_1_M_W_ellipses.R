source("../functions_R/figure_tools.R")

#####################
sims.dirs <- list.dirs("../../simul/fig_1")
#angle = c(-0.2, 1.2)
angle = c(1,2,3) #c(-1.4, -1.3, -1.2 , -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4)
of        <- "fig1_modulopi"
modulo <- pi
#####################


M.factor <- 4
G.factor <- 0.25
#S.factor <- 0.022 #Multilinear
#S.factor <- 0.01 #Multilinear scaled 0_1
#S.factor <- 0.0008 #Wagner
S.factor <- 0.0044

pdfname   <- print(sprintf("../../figures/%s_ellipse.pdf", of))
#Draw ellipses
pdf(pdfname, width=13, height=4.4)
layout(t(1:3))

j <- 1
for (i in angle) {
  j <- j+1
      sims.dirs.w <- list.dirs("../../simul/fig_1/w2")
      sims.dir.w <- sims.dirs.w[j]
      mse.w <- round(mse.ms(diff.ms.df(sims.dir.w, modulo=modulo)), 4)
      sims.dirs.m <- list.dirs("../../simul/fig_1/m_2000")
      sims.dir.m <- sims.dirs.m[j]
      mse.m <- round(mse.ms(diff.ms.df(sims.dir.m, modulo=modulo)), 4)

  sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle__%s.par", i), full.names=TRUE)
  
  # oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor,asp=1, main=sprintf("MSE m = %s MSE w = %s", mse.m, mse.w),xlim=c(-0.5,0.5), ylim=c(-0.5,0.5), all.reps=TRUE, Gell=FALSE) #xlim=c(-3,3), ylim=c(-3,3) all.gen=TRUE all.reps=TRUE
  oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor,asp=1, main=sprintf("MSE m = %s MSE w = %s", mse.m, mse.w),xlim=c(-0.5,0.5), ylim=c(-0.5,0.5), all.reps=TRUE)
  }
dev.off()
print("Ellipses done !")
