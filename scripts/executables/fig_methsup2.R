fig.width  <- 4
fig.height <- 4

col.kappa <- c(orange=0.25, black=0.5, blue=0.75)

f <- function(x, kappa) 1/(1+(1/kappa-1)*exp(-x/kappa/(1-kappa)))

pdf("figures/fig_methsup2.pdf", width=fig.width, height=fig.height)
	par(mar=c(5,4,0.5,0.5))
	
	plot(NULL, xlab="x: Total strength of regulation on gene i", ylab="f(x): Future expression of gene i", xlim=c(-1,1), ylim=c(0,1))
	abline(v=0, lty=2, col="gray")
	
	for (i in seq_along(col.kappa)){
		curve(f(x, col.kappa[i]), add=TRUE, col=names(col.kappa)[i], lwd=2)
		arrows(x0=-0.25, x1=0.25, y0=col.kappa[i]-0.25, y1=col.kappa[i]+0.25, length=0, col=names(col.kappa)[i])
	}
		
	legend("topleft", lty=1, col=names(col.kappa), legend=as.expression(sapply(col.kappa, function(k) bquote(kappa == .(k)))))	
dev.off()

