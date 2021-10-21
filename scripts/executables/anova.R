source("../functions_R/network.R")
source("../functions_R/tools.R")
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(matrixStats)
library(viridis)
#######################
sims.dirs <- c(     
  # "../../simul/wagner/4g/c-1_-1","../../simul/wagner/4g/c1_1","../../simul/wagner/4g/c0_0",
  # "../../simul/wagner/8g/c-1_-1","../../simul/wagner/8g/c1_1","../../simul/wagner/8g/c0_0",
  # "../../simul/wagner/12g/c-1_-1","../../simul/wagner/12g/c1_1","../../simul/wagner/12g/c0_0",
  # "../../simul/wagner/16g/c-1_-1","../../simul/wagner/16g/c1_1","../../simul/wagner/16g/c0_0",
  # "../../simul/wagner/4g_pop500/c0_0","../../simul/wagner/4g_pop500/c1_1","../../simul/wagner/4g_pop500/c-1_-1",
  # "../../simul/wagner/4g_pop10000/c0_0","../../simul/wagner/4g_pop10000/c1_1","../../simul/wagner/4g_pop10000/c-1_-1"
  "../../simul/wagner/4g_cn_ab","../../simul/wagner/4g_cp_ab","../../simul/wagner/4g_fn_ab",
  "../../simul/wagner/4g_fp_ab","../../simul/wagner/4g_no_ab","../../simul/wagner/4g_chn_ab",
  "../../simul/wagner/4g_chp_ab","../../simul/wagner/4g_fhn_ab","../../simul/wagner/4g_fhp_ab"
)
of        <- "controlled_W"
where     <- "wagner"
sims.dir  <- sims.dirs[2:30]
#######################

#pdf(sprintf("../../figures/%s/%s_ANOVA.pdf", where, of), width=6, height=6)
for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:30]
  
  df.topo <- df.topo.raw(sims.dir)
  dfx <- df.topo[,c(9,10,11,12,14,15,16,17,19,20,21,22)]
  dfy <- as.data.frame(df.topo$ang_M)
  df <-as.data.frame(cbind(dfx,dfy))
  names(df)[names(df) == "df.topo$ang_M"] <- "M_angle"
  df[] <- lapply( df, as.numeric)

  fit <- anova(df) # y est la variable numérique et A indique les groupes
  summary(fit)
  
}