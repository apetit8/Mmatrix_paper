source("../functions_R/figure_tools.R")

#####################
sims.dir <- list.dirs("../../simul/Wagner/angle/4locus/angle_0.5_0.5")
sims.dir <- sims.dir[2:9]
of       <- "c0.5_opt0.5"
txtname  <- print(sprintf("../../figures/%s_matrices.txt", of))
#####################


#Write matrices
sink(txtname)
for (dd in sims.dir) {print(dd)
  write.matrix(dd)
}
sink()
print("Matrices done !")