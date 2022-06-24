source("../functions_R/All_functions.R")
#
ifelse(!dir.exists(file.path("../../simul")), dir.create("../../simul"), FALSE)
#General properties of S :
def.e <-  sqrt(1-(0.12)) #eccentricity S
def.s <- 10 #size S
#
################################################################################
#Print paramfiles Figure 1
#####################
values <- list(-pi/8, pi/8, 0, pi/4, 3*pi/8, 0) #Angle values
#####################

ifelse(!dir.exists(file.path("../../simul/fig_1abc")), dir.create("../../simul/fig_1abc"), FALSE)

#This loop creates a parameter template for each angle given.
sims.dirs <- c("../../templates/fig_1abc/m") #multilinear template location
for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "../../templates/", n=2, simplify = TRUE)
  dir <- sprintf("../../simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Multilinear paramfiles fig 1 : done !")


sims.dirs <- c("../../templates/fig_1abc/w") #Wagner template location
for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "../../templates/", n=2, simplify = TRUE)
  dir <- sprintf("../../simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Wagner paramfiles fig 1 : done !")

################################################################################
#Print paramfiles Figure 2
#####################
sims.dirs <- list.dirs("../../templates/fig_1de", recursive = FALSE)
#####################

#This loop creates a parameter template for each angle given.

#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("../../simul/fig_1de")), dir.create("../../simul/fig_1de"), FALSE)

for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "../../templates/", n=2, simplify = TRUE)
  dir <- sprintf("../../simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}

################################################################################
#Print paramfiles Figure 3
#####################
sims.dirs <- list.dirs("../../templates/fig_3", recursive = FALSE)
#####################

#This loop creates a parameter template for each angle given.

#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("../../simul/fig_3")), dir.create("../../simul/fig_3"), FALSE)

for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "../../templates/", n=2, simplify = TRUE)
  dir <- sprintf("../../simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Parameter files fig 3 done !")

################################################################################
#Print paramfiles Figure 4
#####################
sims.dirs <- list.dirs("../../templates/fig_4", recursive = FALSE)
#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("../../simul/fig_4")), dir.create("../../simul/fig_4"), FALSE)


for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "../../templates/", n=2, simplify = TRUE)
  dir <- sprintf("../../simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Parameter files witness done !")

################################################################################
