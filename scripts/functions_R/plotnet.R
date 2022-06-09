
plotnet <- function(W, gene.names=colnames(W), thresh=0, max=1, col.scale.plus=NULL, col.scale.minus=NULL, lwd.arr=2, ...) {
	
	circ.arc    <- function(theta1=0, theta2=2*pi, n=100) 
		{ tt <- seq(theta1, theta2, length.out=n); cbind(cos(tt), sin(tt)) }
	posit.angle <- function(angle)  
		{ angle <- angle %% (2*pi); if (angle > pi/2 && angle <= 3*pi/2) angle <- angle + pi; angle %% (2*pi)}
		
	lg <- ncol(W)
	W[abs(W) < thresh] <- 0
	W[abs(W) > max]    <- max*sign(W[abs(W) > max])
	
	if (is.null(gene.names))
		gene.names <- LETTERS[1:lg]
	if (is.null(col.scale.plus))
		col.scale.plus <- colorRampPalette(c("white","black"))(100)
	if (is.null(col.scale.minus))
		col.scale.minus <- colorRampPalette(c("white","red"))(100)
	
	arr.dist    <- 0.15 # distance between the group name and the arrows
	delta.angle <- 0.25 # angle between two arrows
	self.angle  <- 1.4*pi
	
	xy.genes <- cbind(cos(2*pi/lg*(0:(lg-1))), sin(2*pi/lg*(0:(lg-1))))
	
	par(mar=c(0.1,0.1,2,0.1))
	plot(NULL, xlim=c(-1.2,1.2), ylim=c(-1.2,1.2), axes=FALSE, ann=FALSE, asp=1, ...)
	
	text(xy.genes[,1], xy.genes[,2], parse(text=gene.names))
	
	for (i in 1:lg) {
		for (j in 1:lg) {
		  
			col <- if(W[j,i] > 0) col.scale.plus [round(abs(W[j,i]) * length(col.scale.plus ))]
			       else           col.scale.minus[round(abs(W[j,i]) * length(col.scale.minus))]
			       
			if (i != j) {
				# angle between i and j
				alpha <- atan((xy.genes[j,2]-xy.genes[i,2])/(xy.genes[j,1]-xy.genes[i,1]))
				if (xy.genes[i,1] > xy.genes[j,1]) 
					alpha <- alpha - pi
				
				x.i  <- xy.genes[i,1]+arr.dist*cos(alpha-delta.angle/2)
				y.i  <- xy.genes[i,2]+arr.dist*sin(alpha-delta.angle/2)
				x.j  <- xy.genes[j,1]+arr.dist*cos(pi+alpha+delta.angle/2)
				y.j  <- xy.genes[j,2]+arr.dist*sin(pi+alpha+delta.angle/2)
				
				arrows(x0=x.i,  x1=x.j,  y0=y.i,  y1=y.j,  length=0.1,  col=col,  lwd=lwd.arr)
			} else { # i == j
				# angle of i around the circle
				alpha <- (i-1)*2*pi/lg
				# the regulation circle is opposite to the center
				cc <- circ.arc(alpha-self.angle/2, alpha+self.angle/2)
				cc <- t(t(cc)*arr.dist+(1+0.5*arr.dist)*xy.genes[i,])
				lines(cc[1:(nrow(cc)-1), 1], cc[1:(nrow(cc)-1),2], col=col, lty=1, lwd=lwd.arr)
				arrows(x0=cc[nrow(cc)-1,1], x1=cc[nrow(cc),1], y0=cc[nrow(cc)-1,2], y1=cc[nrow(cc),2], col=col, length=0.1, lwd=lwd.arr)
			}
		}
	}
}
