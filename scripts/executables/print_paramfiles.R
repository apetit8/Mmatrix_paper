source("scripts/functions_R/All_functions.R")
#
ifelse(!dir.exists(file.path("simul")), dir.create("simul"), FALSE)
#General properties of S :
def.e <-  sqrt(1-(0.12)) #eccentricity S
def.s <- 10 #size S
#
################################################################################
#Print paramfiles Figure 1, panels A, B, C
#####################
values <- list(-pi/8, pi/8, 0, pi/4, 3*pi/8, 0) #Angle values
#####################

ifelse(!dir.exists(file.path("simul/fig_2ab")), dir.create("simul/fig_2ab"), FALSE)

#This loop creates a parameter template for each angle given.
sims.dirs <- c("templates/fig_2ab/m") #multilinear template location
for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "templates/", n=2, simplify = TRUE)
  dir <- sprintf("simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e, sel.genes=0)
  }
}
print("Multilinear paramfiles fig 2 : done !")

sims.dirs <- c("templates/fig_2ab/w") #Wagner template location
for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "templates/", n=2, simplify = TRUE)
  dir <- sprintf("simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Wagner paramfiles fig 1 : done !")

sims.dirs <- c("templates/fig_2ab/fkl") #Wagner template location
for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "templates/", n=2, simplify = TRUE)
  dir <- sprintf("simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e, sel.genes=0)
  }
}
print("FKL paramfiles fig 2ab : done !")

################################################################################
#Print paramfiles Figure 1 panel D, E + Figure 2 (PLS)
#####################
sims.dirs <- list.dirs("templates/fig_2cd", recursive = FALSE)
#####################

#This loop creates a parameter template for each angle given.

#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("simul/fig_2cd")), dir.create("simul/fig_2cd"), FALSE)
sel.genes <- 0
i <- 1
for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "templates/", n=2, simplify = TRUE)
  dir <- sprintf("simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e, sel.genes=sel.genes)
  }
  i <- i+1
  if(i==3){sel.genes <- c(3:(2+(6-2)/2))}
  }
print("FKL paramfiles fig 2cd & 3: done !")

################################################################################
#Print paramfiles Figure 3
#####################
sims.dirs <- list.dirs("templates/fig_4", recursive = FALSE)
#####################

#This loop creates a parameter template for each angle given.

#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("simul/fig_4")), dir.create("simul/fig_4"), FALSE)

for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "templates/", n=2, simplify = TRUE)
  dir <- sprintf("simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Parameter files fig 4 done !")

################################################################################
#Print paramfiles Figure 4
#####################
sims.dirs <- list.dirs("templates/fig_5", recursive = FALSE)
#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("simul/fig_5")), dir.create("simul/fig_5"), FALSE)


for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "templates/", n=2, simplify = TRUE)
  dir <- sprintf("simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e)
  }
}
print("Parameter files fig 5 done !")

################################################################################
################################################################################
#Print paramfiles Sypplementary Figure 2
#####################
sims.dirs <- list.dirs("templates/fig_supp4", recursive = FALSE)
#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("simul/fig_supp4")), dir.create("simul/fig_supp4"), FALSE)


for (sims.dir in sims.dirs) {
  param.template = file.path(sims.dir, "template.temp")
  new.dir<- str_split(sims.dir, "templates/", n=2, simplify = TRUE)
  dir <- sprintf("simul/%s", new.dir[2])
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=def.s, eccentricity=def.e, sel.genes=c(3:4))
  }
}
print("Parameter files fig supp 4 done !")

################################################################################
################################################################################
#Print paramfiles Param explo
#####################
sims.dirs <- list.dirs("templates/param_explo", recursive = FALSE)
# sims.dirs <- "templates/param_explo/mutsd"
#####################
#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("simul/param_explo")), dir.create("simul/param_explo"), FALSE)

param.template <- list.files(path = sims.dirs, full.names=TRUE, pattern = "\\.temp$")

for (temp in param.template) {
  dir <-  str_split(str_split(temp, "/", simplify = TRUE)[4], "\\.temp", simplify = TRUE)[1]
  dir <- paste0("simul/param_explo/", dir)
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(temp, param.file, angle=a, size=def.s, eccentricity=def.e, sel.genes=c(3:4))
  }
}

print("Parameter files param explo done !")



################################################################################
################################################################################
#Print paramfiles S explo
#####################
param.template <- "templates/fig_1de/3-grn/template.temp"
#####################

#This loop creates a parameter template for each angle given.
#Angle values :
values <- list(-1.5, 1.5,-1.4, 1.4, -1, 1, -0.5, 0.5, -0.7, -0.2, 0.2, -0.3, 0.3, -0.4, 0.4, 0, -0.6, 0.6, -0.9, 0.9, 0.7, -0.8, 0.8, -1.1, 1.1, -1.2, 1.2,-1.3, 1.3, -1.4, 1.4,-0.1, 0.1)

ifelse(!dir.exists(file.path("simul/param_explo")), dir.create("simul/param_explo"), FALSE)

s.size <- c(0.1,1,10,100,1000) #size S

for (size in s.size) {
  dir <- paste0("simul/param_explo/Ssize_", size)
  dir.create(dir)
  for (a in values) {
    param.file <- file.path(dir, sprintf("simuangle%s.par", round(a, digits=3)))
    param.from.sel.features(param.template, param.file, angle=a, size=size, eccentricity=def.e)
  }
}

print("Parameter files S size explo done !")

################################################################################
################################################################################
#Print paramfiles multi correlations
#####################
template    <- "./templates/fig_6cor/template.par"
simu.path   <- "./simul/fig_6cor"
#####################
reps    <- 100
netsize <- c(4,10,20,30)
part.sel<- 1.

if (!dir.exists(simu.path)) dir.create(simu.path)

par2launch <- NULL

for (ns in netsize) {
  dns <- file.path(simu.path, paste0("net-n", ns))
  if (!dir.exists(dns)) dir.create(dns)
  
  myparns <- read.param(template)
  
  myparns$GENET_NBLOC       <- ns
  myparns$GENET_MUTRATES    <- ns*0.00166666667/0.1 # mutsd * mutrate / nbr genes = 0.00166666667 
  myparns$INIT_ALLELES_FULL <- rep(0, ns*ns)
  myparns$TYPE_ALLELES      <- c(ifelse(diag(ns) == 0, "normal", "immut"))
  myparns$FITNESS_STRENGTH  <- c(rep(0.05, floor(part.sel*ns)), rep(0, (ns-floor(part.sel*ns))))
  
  for (rr in 1:reps) {
    mypar <- myparns
    
    mypar$FITNESS_OPTIMUM     <- rep(0, ns) # runif(ns, 0, 1)
    
    cormat <- tcrossprod(runif(ns, -1, 1)) # Ensures that it is positive definite (but correlations are not uniform)
    mypar$FITNESS_CORRELATION <- c(cormat[upper.tri(cormat)])
    
    myparfile <- file.path(dns, paste0("replicate-", formatC(rr, width=floor(log10(reps))+1, flag="0")))
    write.param(mypar, paste0(myparfile, ".par"))
    par2launch <- append(par2launch, path.expand(myparfile))
  }
}
print("Parameter files fig_6cor done !")