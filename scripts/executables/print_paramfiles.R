source("../functions_R/All_functions.R")
################################################################################
#Print paramfiles Figure 1
#####################
#Properties of S for multilinear model :
def.e <-  0.12 #eccentricity S
def.s <- 10 #size S
values <- list(-pi/8, pi/8, 0, pi/4, 3*pi/8, -1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1) #Angle values
#####################

ifelse(!dir.exists(file.path("../../simul/fig_1")), dir.create("../../simul/fig_1"), FALSE)

#This loop creates a parameter template for each angle given.
sims.dirs <- c("../../templates/fig_1/m") #multilinear template location
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


sims.dirs <- c("../../templates/fig_1/w") #Wagner template location
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
sims.dirs <- c("../../templates/fig_2/1-mult", "../../templates/fig_2/2-full")

#Properties of S for Wagner model :
def.e <- 0.12 #1 #eccentricity S
def.s <- 10 #size S

#####################

#This loop creates a parameter template for each angle given.

#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("../../simul/fig_2")), dir.create("../../simul/fig_2"), FALSE)

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
#
#####################
sims.dirs <- c("../../templates/fig_2/3-round_S")

#Properties of S for Wagner model :
def.e <- 1 #1 #eccentricity S

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
print("Parameter files fig 2 done !")

################################################################################
#Print paramfiles Figure 2
#####################
sims.dirs <- list.dirs("../../templates/fig_3", recursive = FALSE)

#Properties of S for Wagner model :
def.e <- 0.12 #1 #eccentricity S
def.s <- 10 #size S

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
sims.dirs <- list.dirs("../../templates/witness", recursive = FALSE)
#Properties of S for Wagner model :
def.e <- 1 #1 #eccentricity S
def.s <- 10 #size S

#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("../../simul/witness")), dir.create("../../simul/witness"), FALSE)


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
