source("scripts/functions_R/All_functions.R")
#General properties of S :
def.e <-  sqrt(1-(0.12)) #eccentricity S
def.s <- 10 #size S
#
################################################################################
#Print paramfiles Figure Supp
#####################
df.fig1de <- df.data("simul/fig_1de/3-grn", pattern = "simul/fig_1de/", variable="netw", file_size=15000, w_of_6=TRUE, network=TRUE)
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)
#####################
#Lsit of angle to reach
#In all of my full netw simulations, which is closest to 
#Print it in paramfile
#To launch after all others

dir <- "simul/fig_supp6_initW/grn"
ifelse(!dir.exists(file.path(dir)), dir.create(dir), FALSE)
param.template = file.path("templates/fig_1de/3-grn", "template.temp")
for (a in values) {
  #Find closest network
  mvalue <- which(abs(df.fig1de$ang_M - a) == min(abs(df.fig1de$ang_M - a)))
  W <- c(unlist( t(matrix(df.fig1de[mvalue,11:46], ncol=6)), use.names = FALSE))

  #Put S matrix in template file
  print("Parameter files fig_supp6 done !")
  param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
  param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e, initW=W)
}

