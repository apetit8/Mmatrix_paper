source("../functions_R/mutsd.R")

simevolv =("../../../simevolv/bin/Release/Simul_Prog")

sims.dir       = ("../../simul/Multilinear/mutsd")

multi.template = file.path(sims.dir, "multilinear.temp")
wagner.template = file.path(sims.dir, "wagner.temp")

calibrate.wag.on.multi(param.multilin=multi.template, param.wagner=wagner.template)



calibrate.multilin(param.multilin=multi.template, param.wagner=wagner.template)



system(paste(sim.path, "-p", pff, "-o", off, sep=" "))





#dd <- dd[grep(dd, pattern="c0.5*")] # just in case there are directories without simulations
#mutsd.max <- 0.1

# M.fulldat <- as.data.frame(do.call(rbind, lapply(dd, function(ddd) {
# 	param.file <- list.files(path=ddd, pattern="*.par", full.names=TRUE)
# 	param <- read.param(param.file) # from tools.R
# 	mut.sd <- param$GENET_MUTSD
# 	sim.files <- list.files(path=ddd, pattern="simul.*.txt", full.names=TRUE)
# 	Mmat <- lapply(sim.files, function(ff) c(unlist(extract.M.matrix(read.table(ff, header=TRUE)))))
# 	data.frame(mutsd=rep(mut.sd, length(Mmat)), Mmat=do.call(rbind, Mmat))
# })))
# 
# pdf("mutsd-figs.pdf", width=10, height=5)
# 	layout(t(1:2))
# 	plot(M.fulldat$mutsd, M.fulldat$Mmat.1, log="xy", xlab="Mut.sd", ylab="M(1,1)")
# 	plot(log(M.fulldat$mutsd), log(M.fulldat$Mmat.1), xlab="log(Mut.sd)", ylab="log(M(1,1))")
# 	ll <- lm(log(M.fulldat$Mmat.1)[M.fulldat$mutsd <= mutsd.max] ~ log(M.fulldat$mutsd[M.fulldat$mutsd <= mutsd.max]))
# 	abline(ll)
# 	text(x=median(log(M.fulldat$mutsd))+2, y=median(log(M.fulldat$Mmat.1))-1, labels=paste0("y = ", round(coef(ll)[2], digits=2), "x + ", round(coef(ll)[1], digits=2)))
# dev.off()
# 
# 
# print(ll)
# 
# mean(log(M.fulldat$Mmat.1))
# mean(M.fulldat$Mmat.1)
