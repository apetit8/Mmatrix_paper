fig.width  <- 4
fig.height <- 4

col.M <- "red"

pdf("figures/fig_methsup1.pdf", width=fig.width, height=fig.height)
	par(mar=c(2,2,0.5,0.5))
	plot(NULL, xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n", xlab="Trait a", ylab="Trait b", asp=1, mgp=c(1,1,1))
	
	abline(h=0)
	
	lines(ellipse::ellipse(0.8, scale=c(1, 1), t=1), col=col.M, lwd=3)
	text(-0.2, 0.7, "M", col=col.M, cex=1.4)
	arrows(x0=-1, y0=-1, x1=1, y1=1, length=0, col=col.M)
	
	circ <- ellipse::ellipse(0, scale=c(0.2, 0.2), t=1)
	circ <- circ[circ[,1] > 0.2*sqrt(2)/2 & circ[,2] > 0,]
	lines(circ)
	arrows(x0=circ[nrow(circ)-1,1],x1=circ[nrow(circ),1],y0=circ[nrow(circ)-1,2],y1=circ[nrow(circ),2], length=0.1)
	text(0.2,0.1,expression(alpha(M)), pos=4)

dev.off()
