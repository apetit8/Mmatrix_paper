source("../functions_R/figure_tools.R")

#####################
sims.dirs <- c(     
                    # "../../simul/wagner/4g_c1_1_Wneg" 
                    # "../../simul/wagner/4g/c-1_-1","../../simul/wagner/4g/c1_1","../../simul/wagner/4g/c0_0",
                    # "../../simul/wagner/8g/c-1_-1","../../simul/wagner/8g/c1_1","../../simul/wagner/8g/c0_0",
                    # "../../simul/wagner/12g/c-1_-1","../../simul/wagner/12g/c1_1","../../simul/wagner/12g/c0_0",
                    # "../../simul/wagner/16g/c-1_-1","../../simul/wagner/16g/c1_1","../../simul/wagner/16g/c0_0"
                    # "../../simul/wagner/4g_pop500/c0_0","../../simul/wagner/4g_pop500/c1_1","../../simul/wagner/4g_pop500/c-1_-1",
                    # "../../simul/wagner/4g_pop10000/c0_0","../../simul/wagner/4g_pop10000/c1_1","../../simul/wagner/4g_pop10000/c-1_-1"
                    # "../../simul/wagner/4g_cn_ab","../../simul/wagner/4g_cp_ab","../../simul/wagner/4g_fn_ab",
                    # "../../simul/wagner/4g_fp_ab","../../simul/wagner/4g_no_ab","../../simul/wagner/4g_chn_ab",
                    # "../../simul/wagner/4g_chp_ab","../../simul/wagner/4g_fhn_ab","../../simul/wagner/4g_fhp_ab"
                    # "../../simul/fig_2/c0_0_free","../../simul/fig_2/c0_0abba","../../simul/fig_2/c0_0adda", "../../simul/fig_2/round_s"
  "../../simul/fig_2/c0_0_free","../../simul/fig_2/c0_0abba","../../simul/fig_2/c0_0adda", "../../simul/fig_2/c0_0ab"
  )
of        <- "fig2_Ow_modulopi"
what      <- "angle"
modulo <-pi
#####################

sims.dirs.wit <- list.dirs("../../simul/fig_2/witness_0w")
sims.dir.wit <- sims.dirs.wit[2:32]
diff.wit <- diff.ms.df(sims.dir.wit, modulo = modulo)

pdfname   <- print(sprintf("../../figures/%s_MS_plot.pdf", of))
pdf(pdfname, width=9, height=10)
  # layout(t(2:2))
par(mfrow = c(2, 2))
  
  for (i in sims.dirs) {
    sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
    diff <- diff.ms.df(sims.dir, modulo = modulo)
    mse <- mse.ms(diff)
    ks <- ks.test(((diff)^2),((diff.wit)^2))
    wx <- wilcox.test(((diff)^2),((diff.wit)^2))
    sdks <- ""
    sdwx <- ""
    if (ks$p.value < 0.05) {sdks <- "*"}
    if (wx$p.value < 0.05) {sdwx <- "*"}
    pop <- str_split(i, "../../simul/fig_", n=2, simplify = TRUE)
    plot.features.onS(sims.dir, what=what, main=sprintf("%s MSE = %s ks %s wx %s", pop[1,2], round(mse, 4), sdks, sdwx), generation =10000,
                      all.reps=TRUE, xlim =c(-pi/2,pi/2) , ylim=c(-pi/2,pi/2), asp=1, axes=FALSE )
    abline(coef = c(0,1), col="orange")
    axis(side=2, at=c(-pi, -pi/2, -pi/4, 0 ,pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi))
    axis(side=1, at=c(-pi, -pi/2, -pi/4, 0 ,pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi))
  }
    
dev.off()
print("Plot /S done !")