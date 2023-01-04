#!/usr/bin/env Rscript

source("./scripts/functions_R/All_functions.R")

template    <- "./templates/fig_cor/template.par"
simu.path   <- "./simul/fig_cor"
launch.path <- "./launcher/fullcor-launch.sh"
simul.prog  <- file.path(Sys.getenv("HOME"), "Work/Software/simevolv/bin/Release/Simul_Prog")

reps    <- 100
netsize <- c(2,5,10,50)
part.sel<- 1.

################# Creating simulation and parameters files #############
if (!dir.exists(simu.path)) dir.create(simu.path)

par2launch <- NULL

for (ns in netsize) {
	dns <- file.path(simu.path, paste0("net-n", ns))
	if (!dir.exists(dns)) dir.create(dns)
	
	myparns <- read.param(template)
	
	myparns$GENET_NBLOC       <- ns
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

# ################ Creating the launch file ##############################
# 
# launch <- sapply(par2launch, function(pp) {
# 		paste(simul.prog, "-p", paste0(pp, ".par"), "-o", paste0(pp, ".txt"))
# 	})
# 	
# cat(launch, file=launch.path, sep="\n")
