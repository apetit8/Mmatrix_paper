source("../functions_R/mutsd.R")
#################################
simevolv =("../../../simevolv/bin/Release/Simul_Prog")
sims.dir       = ("../../templates/fig_1")
#################################

multi.template = file.path(sims.dir, "m/template.temp")
wagner.template = file.path(sims.dir, "w/template.temp")
calibrate.multilin(param.multilin=multi.template, param.wagner=wagner.template)
