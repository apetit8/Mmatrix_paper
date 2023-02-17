source("scripts/functions_R/All_functions.R")

library(viridis)

#######################
sim.dir      <- "simul/fig_6cor"
#######################

net.size.dir <- list.dirs(sim.dir, full.names=FALSE, recursive=FALSE)
# net.sizes <- as.numeric(gsub("net-n", "", net.size.dir))
net.sizes <- c(10, 30)
net.size.dir <- net.size.dir[order(net.sizes)]
net.sizes    <- net.sizes[order(net.sizes)]

max.dots <- 1000

get.fitcor <- function(parfile, only.sel=TRUE) {
	pp <- read.param(parfile)
	if (only.sel) {
		sel <- which(pp$FITNESS_STRENGTH > 0.01)
	} else {
		sel <- seq_len(pp$GENET_NBLOC)
	}
	cormat <- diag(pp$GENET_NBLOC)
	cormat[upper.tri(cormat)] <- pp$FITNESS_CORRELATION
	cormat <- cormat[sel,sel]
	c(cormat[upper.tri(cormat)])
}

get.mutcor <- function(resfile, parfile, gen=NULL, only.sel=TRUE) {
	pp <- read.param(parfile)
	if (only.sel) {
		sel <- which(pp$FITNESS_STRENGTH > 0.01)
	} else {
		sel <- seq_len(pp$GENET_NBLOC)
	}
	mm <- if (is.null(gen)) extract.M.matrix(read.table(resfile, header=TRUE)) else extract.M.matrix(read.table(resfile, header=TRUE), gen=gen)
	cormat <- cov2cor(mm)[sel,sel]
	c(cormat[upper.tri(cormat)]) # or lower.tri? 
}

get.mutecc <- function(resfile, parfile, gen=NULL, only.sel=TRUE) {
	pp <- read.param(parfile)
	if (only.sel) {
		sel <- which(pp$FITNESS_STRENGTH > 0.01)
	} else {
		sel <- seq_len(pp$GENET_NBLOC)
	}
	mm <- if (is.null(gen)) extract.M.matrix(read.table(resfile, header=TRUE)) else extract.M.matrix(read.table(resfile, header=TRUE), gen=gen)
	
	eccmat <- matrix(1, ncol=length(sel), nrow=length(sel))
	
	for (i1 in seq_along(seq)) {
		for (i2 in seq_along(seq)) {
			if (i1 < i2) {
				ev <- eigen(m[c(sel[i1],sel[i2]),c(sel[i1],sel[i2])])$values
				eccmat[i1,i2] <- sqrt(1-ev[2]/ev[1])
			}
		}
	}
	c(eccmat[upper.tri(eccmat)])
}

analyze.rep <- function(ff) {
	parfile <- paste0(ff, ".par")
	resfile <- paste0(ff, ".txt")
	
	fitcor <- get.fitcor(parfile)
	
	dd <- read.table(resfile, header=TRUE)
	mutcors <- as.data.frame(lapply(dd$Gen, get.mutcor, resfile=resfile, parfile=parfile))
	
	plot(NULL, xlim=range(dd$Gen), ylim=c(-1,1), xlab="Generations", ylab="Correlation")
	abline(h=fitcor, col=seq_along(fitcor), lty=3)
	for (i in 1:nrow(mutcors)) points(dd$Gen, unlist(mutcors[i,]), pch=1, col=i, type="b")
}

analyze.M <- function(resfile, genes=1:2, extract.FUN=extract.M.matrix) {
	dd <- read.table(resfile, header=TRUE)
	mutcovs <- lapply(dd$Gen, function(gg) extract.FUN(dd, gg)[genes,genes])
	
	plot(NULL, xlim=range(dd$Gen), ylim=range(unlist(mutcovs)), xlab="Generations", ylab="(Co)Variances")
	for (i in 1:length(genes)) {
		lines(dd$Gen, sapply(mutcovs, function(m) m[genes[i],genes[i]]), col=i)
		if (i < length(genes))
			for (j in (i+1):length(genes))
				lines(dd$Gen, sapply(mutcovs, function(m) m[genes[i],genes[j]]), col=i, lty=3)
	}
}

pdf("figures/fig_6cor.pdf", width=7, height=3.5*ceiling(length(net.sizes)))
	eccpalette <- col2rgb(plasma(1000))

	layout(matrix(1:length(net.sizes), byrow=TRUE, ncol=2))
	for (nsi in seq_along(net.sizes)) {
		dpath <- file.path(sim.dir, net.size.dir[nsi])
	
		parfiles <- list.files(path=dpath, pattern=".par", full.names=TRUE, recursive = TRUE)
		resfiles <- list.files(path=dpath, pattern=".txt", full.names=TRUE, recursive = TRUE)
		fitcor <- lapply(seq_along(parfiles), function(i) get.fitcor(parfiles[i]))
		mutcor <- lapply(seq_along(parfiles), function(i) get.mutcor(resfiles[i], parfiles[i]))
		mutecc <- lapply(seq_along(parfiles), function(i) get.mutecc(resfiles[i], parfiles[i]))
		smpl <- sample(seq_along(unlist(fitcor)), min(max.dots, length(unlist(fitcor))))
		plot(unlist(fitcor)[smpl], unlist(mutcor)[smpl], main=paste0("Network of size n=", net.sizes[nsi]),
		     xlab="Fitness correlation r(S)", ylab="Mutational correlation r(M)", xlim=c(-1,1), ylim=c(-1,1), col = rgb(red = eccpalette["red",round(1000*unlist(mutecc))], green = eccpalette["green",round(1000*unlist(mutecc))], blue = eccpalette["blue",round(1000*unlist(mutecc))], alpha = 0.2)) 
		abline(lm(unlist(mutcor) ~ unlist(fitcor)))
	}
dev.off()
