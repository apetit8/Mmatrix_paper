source("../functions_R/All_functions.R")
#####################
# Witness
#####################
sims.dirs <- c(
  "../../simul/fig_3/rround"
)

#Localisation of Simevolv program :
sim <- "/home/apetit/simevolv/bin/Release/Simul_Prog"

#Properties of S for Wagner model :
def.e <- 1 #0.12 #1 #eccentricity
def.s <- 10 #size

#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-3*pi/8, -pi/4, -pi/8, 0, pi/8, pi/4, 3*pi/8, pi/2-0.001)

for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  for (a in values) {
    param.file <- file.path(sims.dir, sprintf("templateangle%s.temp", round(a, digits=3)))
    new.dir<- str_split(param.file, "angle", n=2, simplify = TRUE)
    dir.create(sprintf("%s/%s",sims.dir, new.dir[2]))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
    file.rename(from=param.file, to=sprintf("%s/%s/%s",sims.dir, new.dir[2],sprintf("templateangle%s.temp", round(a, digits=3))))
  }
}

#This loop creates a parameter file for each S angle template
sims.dir <- list.dirs(sims.dirs)[2:length(list.dirs(sims.dirs))]
for (i in sims.dir) {
   templist <- list.files(path=i, pattern='template.*\\.temp', full.names = TRUE)
        for (temp in templist) {
          generate.param.list(template=temp, param.list=list(FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                                                             FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)), vec.indx=c(1,2), simevolv=sim ,reps=13, launchfile="launcher.sh")
        }
}
print("Templates done !")

###############################################################################################################################
#Panel 1
###############################################################################################################################
sims.dirs <- c(
  "../../simul/fig_3/rnul"
)

#Localisation of Simevolv program :
sim <- "/home/apetit/simevolv/bin/Release/Simul_Prog"

#Properties of S for Wagner model :
def.e <- 0.12 #1 #eccentricity
def.s <- 10 #size

#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-3*pi/8, -pi/4, -pi/8, 0, pi/8, pi/4, 3*pi/8, pi/2-0.001)


for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  for (a in values) {
    param.file <- file.path(sims.dir, sprintf("templateangle%s.temp", round(a, digits=3)))
    new.dir<- str_split(param.file, "angle", n=2, simplify = TRUE)
    dir.create(sprintf("%s/%s",sims.dir, new.dir[2]))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
    file.rename(from=param.file, to=sprintf("%s/%s/%s",sims.dir, new.dir[2],sprintf("templateangle%s.temp", round(a, digits=3))))
  }
}

#This loop creates a parameter file for each S angle template
sims.dir <- list.dirs(sims.dirs)[2:length(list.dirs(sims.dirs))]
for (i in sims.dir) {
  templist <- list.files(path=i, pattern='template.*\\.temp', full.names = TRUE)
  for (temp in templist) {
    generate.param.list(template=temp, param.list=list(FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                                                       FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)), vec.indx=c(1,2), simevolv=sim ,reps=13, launchfile="launcher.sh")
  }
}
print("Templates done !")


###############################################################################################################################
#Panel 2
###############################################################################################################################
sims.dirs <- c(
  "../../simul/fig_3/rneg"
)

#Localisation of Simevolv program :
sim <- "/home/apetit/simevolv/bin/Release/Simul_Prog"

#Properties of S for Wagner model :
def.e <- 0.12 #1 #eccentricity
def.s <- 10 #size

#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-3*pi/8, -pi/4, -pi/8, 0, pi/8, pi/4, 3*pi/8, pi/2-0.001)


for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  for (a in values) {
    param.file <- file.path(sims.dir, sprintf("templateangle%s.temp", round(a, digits=3)))
    new.dir<- str_split(param.file, "angle", n=2, simplify = TRUE)
    dir.create(sprintf("%s/%s",sims.dir, new.dir[2]))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
    file.rename(from=param.file, to=sprintf("%s/%s/%s",sims.dir, new.dir[2],sprintf("templateangle%s.temp", round(a, digits=3))))
  }
}

#This loop creates a parameter file for each S angle template
sims.dir <- list.dirs(sims.dirs)[2:length(list.dirs(sims.dirs))]
for (i in sims.dir) {
  templist <- list.files(path=i, pattern='template.*\\.temp', full.names = TRUE)
  for (temp in templist) {
    generate.param.list(template=temp, param.list=list(FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                                                       FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)), vec.indx=c(1,2), simevolv=sim ,reps=13, launchfile="launcher.sh")
  }
}
print("Templates done !")
###############################################################################################################################
sims.dirs <- c(
  "../../simul/fig_3/rpos"
)

#Localisation of Simevolv program :
sim <- "/home/apetit/simevolv/bin/Release/Simul_Prog"

#Properties of S for Wagner model :
def.e <- 0.12 #1 #eccentricity
def.s <- 10 #size

#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-3*pi/8, -pi/4, -pi/8, 0, pi/8, pi/4, 3*pi/8, pi/2-0.001)


for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  for (a in values) {
    param.file <- file.path(sims.dir, sprintf("templateangle%s.temp", round(a, digits=3)))
    new.dir<- str_split(param.file, "angle", n=2, simplify = TRUE)
    dir.create(sprintf("%s/%s",sims.dir, new.dir[2]))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
    file.rename(from=param.file, to=sprintf("%s/%s/%s",sims.dir, new.dir[2],sprintf("templateangle%s.temp", round(a, digits=3))))
  }
}

#This loop creates a parameter file for each S angle template
sims.dir <- list.dirs(sims.dirs)[2:length(list.dirs(sims.dirs))]
for (i in sims.dir) {
  templist <- list.files(path=i, pattern='template.*\\.temp', full.names = TRUE)
  for (temp in templist) {
    generate.param.list(template=temp, param.list=list(FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                                                       FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)), vec.indx=c(1,2), simevolv=sim ,reps=13, launchfile="launcher.sh")
  }
}
print("Templates done !")
