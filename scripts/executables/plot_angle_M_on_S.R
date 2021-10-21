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
                    "../../simul/wagner/4g_cn_ab","../../simul/wagner/4g_cp_ab","../../simul/wagner/4g_fn_ab",
                    "../../simul/wagner/4g_fp_ab","../../simul/wagner/4g_no_ab","../../simul/wagner/4g_chn_ab",
                    "../../simul/wagner/4g_chp_ab","../../simul/wagner/4g_fhn_ab","../../simul/wagner/4g_fhp_ab"
                    )
of        <- "controlled_W"
where     <- "wagner"
what      <- "angle"
#####################

pdfname   <- print(sprintf("../../figures/%s/%s_%s_plot_M_on_S.pdf", where, of, what))

pdf(pdfname, width=4.4, height=5)
  layout(t(1:1))
  
  for (i in sims.dirs) {
    sims.dir <- list.dirs(i)[2:30]
    plot.features.onS(sims.dir, what=what, main=i, generation =10000,
                      all.reps=TRUE, xlim =c(-pi/2,pi/2) , ylim=c(-pi/2,pi/2), asp=1, axes=FALSE )
    abline(coef = c(0,1), col="orange")
    axis(side=2, at=c(-pi, -pi/2, -pi/4, 0 ,pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi))
    axis(side=1, at=c(-pi, -pi/2, -pi/4, 0 ,pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi))
  }
    
dev.off()
print("Plot /S done !")