

fit.cor <- 0.8
fit.size <- 10
theta <- c(0,0)

col.trait <- "blue"

fit.var <- fit.size/2*cbind(c(1, fit.cor),c(fit.cor, 1))
fit.fun <- function(x1, x2) exp(-0.5*(t(c(x1,x2) - theta) %*% solve(fit.var) %*% (c(x1,x2) - theta)))

xx <- seq(-2.5,2.5,length.out=101)

fit <- outer(xx, xx, function(x1, x2) mapply(x1, x2, FUN=fit.fun))

cairo_pdf("../../figures/fig_S3A.pdf", width=4.5, height=5)

	par(mar=c(5,5,7,4))

	image(x=xx, y=xx, z=fit, xaxt="n", yaxt="n", xlab=expression("Expression gene "*italic('a')*" "*(P[a])), ylab=expression("Expression gene "*italic('b')*" "*(P[b])), col=hcl.colors(1024, "YlOrRd", rev = TRUE), asp=1)
	title("Correlated fitness function", line=5)
	abline(h=0, v=0, col="white", lwd=2)
	contour(x=xx, y=xx, z=fit, add=TRUE)
	abline(a=0, b=1, lty=2)
	
	ee <- seq(0.1, 0.9, by=0.2)
	axis(1, at=log(ee/(1-ee)), label=as.character(ee))
	axis(2, at=log(ee/(1-ee)), label=as.character(ee))
	
	axis(3, col=col.trait, col.axis=col.trait)
	mtext(expression("Phenotypic trait "*Z[1]), 3, line=2.5, col=col.trait)
	axis(4, col=col.trait, col.axis=col.trait)
	mtext(expression("Phenotypic trait "*Z[2]), 4, line=2.5, col=col.trait)

dev.off()


makeTransparent<-function(someColor, alpha=70)
{ # from https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

cairo_pdf("../../figures/fig_S3B.pdf", width=5, height=5)
	par(mgp=c(1,1,0))
	contour(x=xx, y=xx, z=fit, asp=1, xlim=c(-1,1), ylim=c(-1.5,1.5), xaxt="n", yaxt="n", col="darkgray", xlab=expression("Expression gene "*italic(a)), ylab=expression("Expression gene "*italic(b)))
	
	title("Competition between genotypes\nwith different mutational properties")
	
	polygon(ellipse::ellipse(cbind(c(1,0),c(0,1)), t=1), border=NA, col=makeTransparent("blue"))
	polygon(ellipse::ellipse(cbind(c(1,0.7),c(0.7,1)), t=1), border=NA, col=makeTransparent("red"))
	
	text(-1.8, 1.4, "Genotype B", col="blue", cex=1, pos=4)
	arrows(-1, 1.2, -0.02, 0.02, col="blue", lwd=2)
	text(1.8, 1.2, "Offspring of\ngenotype B", col="blue", cex=1, pos=2)
	arrows(0.5, 1.2, -0.2, 0.7, col="blue", lwd=2)
	
	text(1.8, -1.4, "Genotype R", col="red", cex=1, pos=2)
	arrows(1.3, -1.2, 0.02, -0.02, col="red", lwd=2)
	text(-1.8, -1.4, "Offspring of\ngenotype R", col="red", cex=1, pos=4)
	arrows(-1.2, -1.1, -0.8, -0.7, col="red", lwd=2)

dev.off()
