source("scripts/functions_R/All_functions.R")

fig.width  <- 5
fig.height <- 5

library(viridis)
col.net <- magma(4)[-4]

yshift <- 0.4

sp <- 36000


w1 <- rbind(c(0,1,-1,0,0,0),c(-1,0,1,0,0,1),c(-2,-0.1,0,0.5,0,0.5),c(0,1,0,0,1,1),c(0,-1,0,0,0,-1),c(0,1,1,-1,0,0))
w2 <- rbind(c(0,1,-1,0,0,0),c(-1,0,1,0,0,1),c(-1,-0.88,0,0.5,0,1),c(0,-1,0,0,1,1),c(0,-1,0,0,0,-1),c(0,1,1,-1,0,0))
w3 <- rbind(c(0,1,-1,0,0,0),c(-1,0,1,0,0,1),c(-1,-0.97,0,0.5,0,1),c(0,-1,0,0,1,1),c(0,-1,0,0,0,-1),c(0,1,1,-1,0,0))

meanvar <- function(w) {
	pp <- pheno.from.W(w, full=FALSE)
	list(mean=pp$sumx/4, var=pp$sumx2/4 - (pp$sumx/4)^2)
}

sumV1 <- sum(meanvar(w1)$var)
sumV2 <- sum(meanvar(w2)$var)
sumV3 <- sum(meanvar(w3)$var)

plotl <- function(w, col.genes="black", xlab="Network time steps", ylab=expression("Gene expression "*P[i]), ...) {
	col.genes <- rep(col.genes, length.out=ncol(w))
	pp <- pheno.from.W(w, full=TRUE)
	plot(NULL, xlim=c(0,ncol(pp$full)-1), ylim=c(0,1), xlab=xlab, ylab=ylab, ...)
	polygon(c(21,24,24,21), c(0,0,1,1), col="lightgray")
	for(i in 1:ncol(w))
		lines(x=0:(ncol(pp$full)-1), pp$full[i,], col=col.genes[i])
}



pdf("figures/fig_methsup3.pdf", width=fig.width, height=fig.height)
	layout(rbind(1:3,rep(4,3)))
	
	par(mar=c(4,1,0.5,0.5), oma=c(0,4,0,0))
	plotl(w1, col.genes=col.net[1], xpd=NA, mgp=c(2,1,0))
	plotl(w2, col.genes=col.net[2], ylab="", yaxt="n", mgp=c(2,1,0)); axis(2, labels=FALSE)
	plotl(w3, col.genes=col.net[3], ylab="", yaxt="n", mgp=c(2,1,0)); axis(2, labels=FALSE)
	
	curve(exp(-sp*x), xlim=c(0, 0.001), ylim=c(0,1+yshift), xlab=expression("Network unstability, "*Sigma[i]*V[i]), ylab=expression("Fitness penalty "*w[stab]), yaxt="n", xpd=NA)
	axis(2, at=c(0,1/2,1))
	abline(h=1, col="darkgray")
	arrows(x0=c(sumV1, sumV2, sumV3), y0=1+yshift, y1=exp(-sp*c(sumV1, sumV2, sumV3))+0.1, length=0.1, lwd=4, col=col.net)
	points(c(sumV1, sumV2, sumV3), exp(-sp*c(sumV1, sumV2, sumV3)), col=col.net, pch=16, cex=3)
dev.off()

