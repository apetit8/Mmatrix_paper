source("../functions_R/mutsd.R")
#################################
simevolv =("../../../simevolv/bin/Release/Simul_Prog")
sims.dir       = ("../../simul/fig_1")
#################################
i <- c(1,2,3)

for (j in i) {
  
  multi.template = file.path(sims.dir, sprintf("m/template%s.temp", j))
  wagner.template = file.path(sims.dir, sprintf("w/template%s.temp", j))
  calibrate.multilin(param.multilin=multi.template, param.wagner=wagner.template)
}