source("../functions_R/figure_tools.R")

#####################
sims.dirs <- list.dirs("../../simul/cercle/multi/c0_0")
sims.dir  <- sims.dirs[2:30]
of        <- "mut_c0_0"
where     <- "cercle"
what      <- "angle"
#####################


pdfname   <- print(sprintf("../../figures/%s/%s_%s_featplot.pdf", where, of, what))

pdf(pdfname, width=4.4, height=5)
layout(t(1:1))
  plot.features.onS(sims.dir, what=what, main=sims.dirs[1], generation =5000,
                    all.reps=TRUE, xlim =c(-pi/2,pi/2) , ylim=c(-pi/2,pi/2), asp=1, axes=FALSE )
  abline(coef = c(0,1), col="orange")
  axis(side=2, at=c(-pi, -pi/2, -pi/4, 0 ,pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi))
  axis(side=1, at=c(-pi, -pi/2, -pi/4, 0 ,pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi))
  
dev.off()
print("Plot /S done !")




#xlim =c(-pi/2,pi/2) , ylim=c(-pi/2,pi/2), asp=1
# #source("../script_R/tools.R")
# def.e <- 0.3
# def.s <- 0.14
# def.a <- pi/4
# #
# layout(matrix(1:15, ncol=5, byrow=TRUE))
# par(mar=c(0.1,0.1,2,0.1))
# 
# value<-c(-pi+0.1, -2.5, -pi/2, -1.3, -0.5, -1, 0, 1, 0.5, 1.3, pi/2, 2.5, pi-0.1)
# for (a in value) {
#   mm <- matrix2.from.features(angle=modulopi(a), size=def.s, eccentricity=def.e)
#   print(mm)
#   plot(NULL, xlim=c(-5,5), ylim=c(-5,5), asp=1, main=paste0("a=",
#                                                             round(a, digits=2)), axes=FALSE)
#   lines(ellipse::ellipse(mm))
#   lines(ellipse.axes(mm))
# }

# for (e in seq(0.001, 0.999, length.out=5)) {
#   mm <- matrix2.from.features(angle=def.a, size=def.s, eccentricity=e)
#   plot(NULL, xlim=c(-5,5), ylim=c(-5,5), asp=1, main=paste0("e=",
#                                                             round(e, digits=2)), axes=FALSE)
#   lines(ellipse::ellipse(mm))
#   lines(ellipse.axes(mm))
# }
# 
# for (s in seq(0.1, 10, length.out=5)) {
#   mm <- matrix2.from.features(angle=def.a, size=s, eccentricity=def.e)
#   plot(NULL, xlim=c(-5,5), ylim=c(-5,5), asp=1, main=paste0("s=",
#                                                             round(s, digits=2)), axes=FALSE)
#   lines(ellipse::ellipse(mm))
#   lines(ellipse.axes(mm))
# }
# 
# 



