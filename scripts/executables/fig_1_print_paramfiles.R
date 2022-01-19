source("../functions_R/All_functions.R")

#####################
#Properties of S for multilinear model :
def.e <-  0.12 #eccentricity
def.s <- 10 #size
values <- list(-pi/8, 0, pi/4, -1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1) #Angle values
#####################

#This loop creates a parameter template for each angle given.
sims.dirs <- c("../../simul/fig_1/m") #multilinear template location
for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  for (a in values) {
    param.file <- file.path(sims.dir, sprintf("templateangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Multilinear paramfiles : done !")


sims.dirs <- c("../../simul/fig_1/w") #Wagner template location
for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  for (a in values) {
    param.file <- file.path(sims.dir, sprintf("templateangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Wagner paramfiles done !")
