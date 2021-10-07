source("../functions_R/tools.R")

#####################
sims.dirs <- c("../../simul/wagner/12g/c-1_-1","../../simul/wagner/12g/c1_1","../../simul/wagner/12g/c0_0")

#Properties of S for Wagner model :
def.e <- 0.12 
def.s <- 0.1
def.a <- pi/4

#####################

#This loop creates a parameter template for each angle given.

#Angle values :
values <- list(-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  for (a in values) {
    param.file <- file.path(sims.dir, sprintf("templateangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Templates done !")