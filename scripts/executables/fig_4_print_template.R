source("../functions_R/tools.R")
######################
list.dirs <- c(
  "../../simul/fig_4/init2"
)
#Properties of S for Wagner model :
def.e <- 0.12 #1 #eccentricity
def.s <- 10 #size
#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-3*pi/8, -pi/8, pi/8, 3*pi/8, 0, pi/4, -pi/4)

sims.dirs <- list.dirs(list.dirs)[2:length(list.dirs(list.dirs))]

for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  for (a in values) {
    param.file <- file.path(sims.dir, sprintf("templateangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
