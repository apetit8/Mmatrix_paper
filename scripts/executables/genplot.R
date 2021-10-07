source("../functions_R/figure_tools.R")

#####################
sims.dirs <- list.dirs("../../simul/Multilinear/jones_test/opt5")
sims.dir  <- sims.dirs[1]
of        <- "opt5"
where     <- "multi/jones"
pdfname   <- print(sprintf("../../figures/%s/%s_genplot.pdf", where, of))
#####################


#Plot matrices features 
pdf(pdfname, width=20, height=5)
layout(t(1:4))
for (dd in sims.dir) {
  print(dd)
  plot.features.time(dd, what = "angle", main=dd)
  plot.features.time(dd, what = "size", main=dd)
  plot.features.time(dd, what = "eccentricity", main=dd)
  plot.features.time(dd, what = "cor", main=dd)
}
dev.off()
print("Plots done !")

features.extract.matrix(sims.dir)
