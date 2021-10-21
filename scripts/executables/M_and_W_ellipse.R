source("../functions_R/figure_tools.R")

#####################
sims.dirs <- list.dirs("../../simul/cercle/test")
#angle = c(-0.2, 1.2)
angle = c(0.1) #c(-1.4, -1.3, -1.2 , -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4)
of        <- "cercle_ex"
where     <- "cercle"
#####################


M.factor <- 2
G.factor <- 0.25
#S.factor <- 0.022 #Multilinear
#S.factor <- 0.01 #Multilinear scaled 0_1
#S.factor <- 0.0008 #Wagner
S.factor <- 0.0022

#pdfname   <- print(sprintf("../../figures/%s/%s_ellipse.pdf", where, of))
#Draw ellipses
#pdf(pdfname, width=7, height=7)
layout(t(1:1))


for (i in angle) {
  
  
  sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle_%s.par", i), full.names=TRUE)
  
  oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, main=i ,
                     xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), all.reps=FALSE, Gell=FALSE)#xlim=c(-3,3), ylim=c(-3,3) all.gen=TRUE all.reps=TRUE
}
#dev.off()
print("Ellipses done !")
