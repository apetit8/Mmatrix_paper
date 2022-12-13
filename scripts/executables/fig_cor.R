source("../functions_R/All_functions.R")

sim.dir      <- "../../simul/fig_cor"
net.size.dir <- list.dirs(sim.dir, full.names=FALSE, recursive=FALSE)

net.sizes <- as.numeric(gsub("net-n", "", net.size.dir))

net.size.dir <- net.size.dir[order(net.sizes)]
net.sizes    <- net.sizes[order(net.sizes)]

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

pdf("../../figures/fig_cor.pdf", width=10, height=10)
	layout(matrix(1:length(net.sizes), byrow=TRUE, ncol=2))
	for (nsi in seq_along(net.sizes)) {
		dpath <- file.path(sim.dir, net.size.dir[nsi])
	
		parfiles <- list.files(path=dpath, pattern=".par", full.names=TRUE)
		resfiles <- list.files(path=dpath, pattern=".txt", full.names=TRUE)
		fitcor <- lapply(seq_along(parfiles), function(i) get.fitcor(parfiles[i]))
		mutcor <- lapply(seq_along(parfiles), function(i) get.mutcor(resfiles[i], parfiles[i]))
		plot(unlist(fitcor), unlist(mutcor), main=paste0("Network of size n=", net.sizes[nsi]),
		     xlab="Fitness correlation", ylab="Mutational correlation", xlim=c(-1,1), ylim=c(-1,1)) 
		abline(lm(unlist(mutcor) ~ unlist(fitcor)))
	}
dev.off()



