source("../functions_R/All_functions.R")
#####################
sims.dirs <- c(
  "../../templates/fig_3_bis/no_ab","../../templates/fig_3_bis/neg_ab","../../templates/fig_3_bis/pos_ab"
)

#Localisation of Simevolv program :
sim <- "/home/apetit/simevolv/bin/Release/Simul_Prog"
#####################

#Second evolutionary phase : correlational selection

#Properties of S:
def.e <- 0.12 #eccentricity
def.s <- 10 #size

##This loop creates a parameter template for each angle given.
values <- list(-3*pi/8, -pi/4, -pi/8, 0, pi/8, pi/4, 3*pi/8, pi/2-0.001)
for (i in sims.dirs) {
  param.template = file.path(i, "nextpar.temp")
  new.dir<- str_split(i, "../../templates/", n=2, simplify = TRUE)
  dir <- sprintf("../../simul/%s", new.dir[2])
  dir.create(dir)
    for (a in values) {
      param.file <- file.path(dir, sprintf("nextpar%s.temp", round(a, digits=1)))
      param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
      file.rename(from=param.file, to=sprintf("%s/%s", dir, sprintf("nextpar%s.par", round(a, digits=1))))
    }
}

###############################
#Fisrt evolutionary phase : getting to the optimum

#Properties of S:
def.e <- 1 #1 #eccentricity

for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "../../templates/", n=2, simplify = TRUE)
  dir <- sprintf("../../simul/%s", new.dir[2])
  nextpar <- list.files(path=dir, pattern='nextpar.*\\.par', full.names = FALSE)
  param.file <- file.path(dir, sprintf("templateangle%s.temp", round(a, digits=3)))
  param.from.sel.features(param.template, param.file, angle=0, size=def.s, eccentricity=def.e)
  file.rename(from=param.file, to=sprintf("%s/%s", dir, sprintf("templateangle%s.temp", 0)))

  ##This loop creates a parameter file for each optimums
  templist <- list.files(path=dir, pattern='templateangle.*\\.temp', full.names = TRUE)
  for (temp in templist) {
    generate.param.list(template=temp, param.list=list(FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                                                       FITNESS_OPTIMUM=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                                                       FILE_NEXTPAR=nextpar), vec.indx=c(1,2,1),
                        simevolv=sim, reps=13, launchfile="launcher.sh")
  }
}