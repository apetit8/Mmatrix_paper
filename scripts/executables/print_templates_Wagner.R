source("../functions_R/tools.R")

#####################
sims.dirs <- c(
  # "../../simul/wagner/4g/c-1_-1","../../simul/wagner/4g/c1_1","../../simul/wagner/4g/c0_0",
  # "../../simul/wagner/8g/c-1_-1","../../simul/wagner/8g/c1_1","../../simul/wagner/8g/c0_0",
  # "../../simul/wagner/12g/c-1_-1","../../simul/wagner/12g/c1_1","../../simul/wagner/12g/c0_0",
  # "../../simul/wagner/16g/c-1_-1","../../simul/wagner/16g/c1_1","../../simul/wagner/16g/c0_0",
  # "../../simul/wagner/4g_pop500/c0_0", "../../simul/wagner/4g_pop500/c1_1", "../../simul/wagner/4g_pop500/c-1_-1",
  # "../../simul/wagner/4g_pop10000/c0_0", "../../simul/wagner/4g_pop10000/c1_1", "../../simul/wagner/4g_pop10000/c-1_-1",
  # "../../simul/wagner/4g_c1_1_Wneg",
  # "../../simul/wagner/4g_cn_ab","../../simul/wagner/4g_cp_ab","../../simul/wagner/4g_fn_ab",
  # "../../simul/wagner/4g_fp_ab","../../simul/wagner/4g_no_ab",
  # "../../simul/wagner/4g_chn_ab","../../simul/wagner/4g_chp_ab","../../simul/wagner/4g_fhn_ab",
  # "../../simul/wagner/4g_fhp_ab", 
  # "../../simul/fig_2/round_s",
  # "../../simul/fig_2/c0_0ab","../../simul/fig_2/c0_0adda","../../simul/fig_2/c0_0abba",
  # "../../simul/fig_2/round_s3"
  # "../../simul/fig_3/se1_ab_neg","../../simul/fig_3/se_1","../../simul/fig_3/se_-1", "../../simul/fig_3/se_0",
  # "../../simul/witness_0w"
                )

#Properties of S for Wagner model :
def.e <- 1  #0.12 #1 #eccentricity
def.s <- 10 #size

#####################

#This loop creates a parameter template for each angle given.

#Angle values :
# values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)
values <- list(-1.5, 1.5)

for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  for (a in values) {
    param.file <- file.path(sims.dir, sprintf("templateangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Templates done !")