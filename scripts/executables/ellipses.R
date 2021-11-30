source("../functions_R/figure_tools.R")

#####################
sims.dirs <- list.dirs(
  "../../simul/fig_3/round_map", recursive = FALSE
)
sims.dir  <- sims.dirs
of        <- "round_s"
what      <- "angle" #"angle" or "ecc"
where     <- "wagner"
#####################

pdfname   <- print(sprintf("../../figures/%s/%s_ellipse.pdf", where, of))
#pdfname   <- print(sprintf("test.pdf", where, of))

M.factor <- 20
#G.factor <- 0.6
G.factor <- 1
S.factor <- 0.02 #Multilinear
#S.factor <- 0.01 #Multilinear scaled 0_1
#S.factor <- 0.0008 #Wagner
#S.factor <- 0.01

#Draw ellipses
pdf(pdfname, width=11, height=5)
  layout(t(1:3))

  for (dd in sims.dir) {
    print(dd)
    # plot.ellipse.dir(dd, S.factor=S.factor, M.factor=M.factor, G.factor=G.factor, main=dd ,xlim=c(-2,0), ylim=c(-2,0), Gell = FALSE)#xlim=c(-3,3), ylim=c(-3,3) all.gen=TRUE all.reps=TRUE
    # plot.ellipse.dir(dd, S.factor=S.factor, M.factor=M.factor, G.factor=G.factor, Gell=FALSE, main=dd ,xlim=c(-3,3), ylim=c(-3,3), all.reps=TRUE)
    # plot.ellipse.dir(dd, S.factor=S.factor, M.factor=M.factor, G.factor=G.factor, Gell=FALSE, main=dd ,xlim=c(-3,3), ylim=c(-3,3), all.reps=FALSE, pattern_simu = "FITNESS_OPTIMUM--1-FITNESS_OPTIMUM2--1.*.txt$" )#xlim=c(-3,3), ylim=c(-3,3) all.gen=TRUE all.reps=TRUE
     plot.ellipse.dir(dd, S.factor=S.factor, M.factor=M.factor, G.factor=G.factor, Gell=FALSE, main=dd ,xlim=c(-3,3), ylim=c(-3,3), all.reps=FALSE, pattern_simu = "FITNESS_OPTIMUM-1-FITNESS_OPTIMUM2-1.*.txt$" )
    # plot.ellipse.dir(dd, S.factor=S.factor, M.factor=M.factor, G.factor=G.factor, main=dd ,xlim=c(-3,3), ylim=c(-3,3), all.gen=TRUE)#xlim=c(-3,3), ylim=c(-3,3) all.gen=TRUE all.reps=TRUE
        }
dev.off()
print("Ellipses done !")

