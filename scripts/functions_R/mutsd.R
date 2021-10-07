library(parallel)

source("../functions_R/tools.R")
sim.path <- "../../../simevolv/bin/Release/Simul_Prog"
can.tests <- 10

M.from.mutsd <- function(param.multilin, mutsd, gen=100) {
	pp <- read.param(param.multilin)
	#stopifnot(pp$TYPE_ARCHI == "multilinear")

	ff <- tempfile()
	pff <- paste0(ff, ".par")
	off <- paste0(ff, ".out")
		
	pp$SIMUL_GENER     <- gen
	pp$SIMUL_OUTPUT    <- gen
	pp$GENET_MUTSD     <- mutsd
	pp$OUT_CANAL_TESTS <- can.tests
	pp$OUT_CANAL_MUTSD <- mutsd
	
	write.param(pp, parfile=pff)
		
	system(paste(sim.path, "-p", pff, "-o", off, sep=" "))
		
	res <- read.table(off, header=TRUE)
	extract.M.matrix(res)
}

M.from.wagner <- function(param.wagner, gen=100) {
	pp <- read.param(param.wagner)
	stopifnot(pp$TYPE_ARCHI %in% c("m2","wagner"))
	
	ff <- tempfile()
	pff <- paste0(ff, ".par")
	off <- paste0(ff, ".out")
		
	pp$SIMUL_GENER     <- gen
	pp$SIMUL_OUTPUT    <- gen
	pp$OUT_CANAL_TESTS <- can.tests
	
	write.param(pp, parfile=pff)
		
	system(paste(sim.path, "-p", pff, "-o", off, sep=" "))
		
	res <- read.table(off, header=TRUE)
	extract.M.matrix(res)	
}

M.from.multi <- function(param.multilin, gen=100) {
  pp <- read.param(param.multilin)
  stopifnot(pp$TYPE_ARCHI %in% c("multilinear"))
  
  ff <- tempfile()
  pff <- paste0(ff, ".par")
  off <- paste0(ff, ".out")
  
  pp$SIMUL_GENER     <- gen
  pp$SIMUL_OUTPUT    <- gen
  pp$OUT_CANAL_TESTS <- can.tests
  
  write.param(pp, parfile=pff)
  
  system(paste(sim.path, "-p", pff, "-o", off, sep=" "))
  
  res <- read.table(off, header=TRUE)
  extract.M.matrix(res)	
}

mutsd.predict <- function(param.multilin, Vm, reps=10, gen=100, mlim=c(0.0001, 1), traits=NA, mc.cores=detectCores()-1) {
	if (is.matrix(Vm)) Vm     <- diag(Vm)
	if (is.na(traits[1])) traits <- seq_along(Vm)
	
	mm <- exp(seq(log(mlim[1]), log(mlim[2]), length.out=reps))
	
	vm <- unlist(mclapply(mm, function(msd) {
		sum(diag(M.from.mutsd(param.multilin, msd, gen=gen))[traits])
	}, mc.cores=mc.cores))
		
	log.mutsd <- log(mm)
	log.Vm    <- log(vm)
	
	bb <- which.max(vm[vm<sum(Vm)])
	if (bb==1) bb <- 2
	if (bb==reps) bb <- reps-1
	log.mutsd <- log.mutsd[(bb-1):(bb+1)]
	log.Vm    <- log.Vm[(bb-1):(bb+1)]
	ll <- lm(log.mutsd ~ log.Vm)
	exp(predict(ll, list(log.Vm=log(sum(Vm)))))
}
	
calibrate.multilin <- function(param.multilin, param.wagner, gen=100, traits=1:2, ...) {
	M.wagner <- M.from.wagner(param.wagner, gen=gen)
	
	Vm <- M.wagner[traits,traits]
	mst <- mutsd.predict(param.multilin, Vm, gen=gen, traits=traits, ...)
		
	pp <- read.param(param.multilin)
	pp$GENET_MUTSD     <- mst
	pp$OUT_CANAL_MUTSD <- mst
	write.param(pp, parfile=param.multilin)
}

calibrate.wag.on.multi <- function(param.multilin, param.wagner, gen=100, traits=1:2, ...) {
  M.multi <- M.from.multi(param.multilin, gen=gen)
  
  Vm <- M.multi[traits,traits]
  mst <- mutsd.predict(param.wagner, Vm, gen=gen, traits=traits, ...)
  
  pp <- read.param(param.wagner)
  pp$GENET_MUTSD     <- mst
  pp$OUT_CANAL_MUTSD <- mst
  write.param(pp, parfile=param.wagner)
}
