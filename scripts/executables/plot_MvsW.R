source("../functions_R/figure_tools.R")

#####################
sims.dirs <- list.dirs("../../simul/Multilinear/eccentricity/l4_0.5_0.2_a0")
sims.dir1 <- sims.dirs[2:11]
sims.dirs <- list.dirs("../../simul/Wagner/eccentricity/l4_0.5_0.2_a0")
sims.dir2 <- sims.dirs[2:11]
of        <- "l4_c0.5_opt0.2_a0"
what      <- "eccentricity"
#####################


pdfname   <- print(sprintf("../../figures/%s/%s_%s_feat_MvsW", what, of, what))

pdf(pdfname, width=7, height=7)
layout(t(1:1))
plot.features.onS.MvsW(sims.dir1, sims.dir2, what=what, main=sims.dirs[1], generation =5000)
abline(coef = c(0,1), col="orange")
dev.off()
print("Plot MvsW done !")

