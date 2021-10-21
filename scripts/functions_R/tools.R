suppressPackageStartupMessages(library(abind))


# Applies the function FUN (typically, mean or var) to each element of the list of matrix/data.frames 
replicate.apply <- function(tt, FUN=mean, ...) {
	# better checking if the data base is consistent
	stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
	arr <- do.call(abind, c(tt, list(along=3))) #list(along ?
	ans <- apply(arr, 1:2, FUN, ...)
	ans
}

# This is faster than replicate.apply(tt, mean)
replicate.mean <- function(tt) {
	stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
	arr <- do.call(abind, c(tt, list(along=3)))
	rowMeans(arr, dims=2)
}

# This is faster than replicate.apply(tt, var)
replicate.var <- function(tt) {
	# This rowVars function comes from https://stat.ethz.ch/pipermail/r-help/2006-April/103001.html
	# Author: David Brahm
	.rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE, twopass=FALSE) {
		if (SumSquares) return(rowSums(x^2, na.rm, dims))
		N <- rowSums(!is.na(x), FALSE, dims)
		Nm1 <- if (unbiased) N-1 else N
		if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
		(rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
	}
	stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
	arr <- do.call(abind, c(tt, list(along=3)))
	.rowVars(arr, dims=2)
}

extract.P.mean <- function(tt, what="MPhen",  gen=tt[nrow(tt), "Gen"]) {
	ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
	if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
	ex
}

extract.matrix <- function(tt, what="CovPhen", gen=tt[nrow(tt), "Gen"]) {
	ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
	if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
	if (sqrt(length(ex)) %% 1 != 0) stop("No way to make a square matrix out of ", length(ex), " elements.")
	matrix(ex, ncol=sqrt(length(ex)))
}

extract.P.matrix <- function(tt, gen=tt[nrow(tt),"Gen"]) {
	extract.matrix(tt, "CovPhen", gen=gen)
}

extract.M.matrix <- function(tt, gen=tt[nrow(tt), "Gen"]) {
	extract.matrix(tt, "MGenCov", gen=gen)
}

extract.W.matrix <- function(tt, gen=tt[nrow(tt), "Gen"]) {
	t(extract.matrix(tt, "MeanAll", gen=gen))
}

##
#function to read a parameter file
read.param <- function(parfile) {
	ff <- readLines(parfile)
	ss <- strsplit(ff, split="\\s+")
	ans <- sapply(ss, function(x) {cx <- suppressWarnings(as.numeric(x[-1])); if(any(is.na(cx))) x[-1] else cx})
	param.names <- sapply(ss, "[", 1)
	if (any(duplicated(param.names))) warning("Duplicated param names in ", parfile)
	names(ans) <- param.names
	ans
}

##
#function to write in a parameter file
write.param <- function(parlist, parfile) {
	unlink(parfile)
	for (nn in names(parlist)) 
		cat(nn, "\t", parlist[[nn]], "\n", file=parfile, append=TRUE)
}

extract.nbphen <- function(parfile) {
	rp <- read.param(parfile)
	if (rp$TYPE_ARCHI %in% c("additive", "multilinear")) {
		np <- if("GENET_NBPHEN" %in% names(rp)) rp$GENET_NBPHEN else 1
	} else if (rp$TYPE_ARCHI %in% c("wagner", "siegal", "m2")) {
		np <- rp$GENET_NBLOC
	}
	np
}

extract.theta <- function(parfile) {
	rp <- read.param(parfile)
	np <- extract.nbphen(parfile)
	th <- rp$FITNESS_OPTIMUM
	if (length(th) == 1) return(rep(th, np))
	if (length(th) == np) return(th)
	stop("Param file: ", parfile, ", The number of optima (", length(th), ") does not match the number of traits (", np, ").")
}

extract.S.matrix <- function(parfile) {
	.cor2cov <- function(cc, sd) t(cc*sd)*sd

	rp <- read.param(parfile)
	np <- extract.nbphen(parfile)
	fs <- rp$FITNESS_STRENGTH
	if (length(fs)==1)
		fs <- rep(fs, np)
	if (length(fs) != np)
		stop("Param file: ", parfile, ", The number of fitness strengths (", length(fs), ") does not match the number of traits (", np, ").")
	if (rp$FITNESS_TYPE == "gaussian") {
		ans <- diag(1/(2*fs))
	} else if (rp$FITNESS_TYPE == "multivar_gaussian") {
		vv  <- 1/(2*fs)
		rr <- rp$FITNESS_CORRELATION
		if (length(rr) != np*(np-1)/2)
			stop("Param file: ", parfile, ", The number of fitness correlations (", length(rr), ") does not match the number of traits (", np, ").")
		ansr <- diag(np)
		ansr[upper.tri(ansr)] <- rr
		ansr <- as.matrix(Matrix::forceSymmetric(ansr))
		ans <- .cor2cov(ansr, sqrt(vv))
	} else {
		stop("Asking the S matrix for non-gaussian selection is meaningless.")
	}
	ans
}

##
#Take properties of S and convert it in S parameters for Simevolv 
param.S.matrix <- function(S.mat, n.genes=ncol(S.mat)) {
	stopifnot(is.matrix(S.mat), ncol(S.mat) == nrow(S.mat))
	stopifnot(n.genes > 0)
	if (n.genes < ncol(S.mat))
		S.mat <- S.mat[1:n.genes, 1:n.genes]
	if (n.genes > ncol(S.mat)) {
		rec <- matrix(0, ncol=n.genes-ncol(S.mat), nrow=nrow(S.mat))
		dg <- matrix(0, ncol=n.genes-ncol(S.mat), nrow=n.genes-nrow(S.mat))
		diag(dg) <- Inf
		S.mat <- rbind(cbind(S.mat, rec), cbind(t(rec), dg))
	}
	list(
		FITNESS_STRENGTH = 1/(2*diag(S.mat)),
		FITNESS_CORRELATION = cov2cor(S.mat)[upper.tri(S.mat)])
}

######### Matrix features

modulopi <- function(angle) {
	# returns the angle of an ellipse (between -pi/2 and pi/2)
	ans <- angle %% pi
	ifelse(ans > pi/2, ans-pi, ans)
}

matrix.features <- function(M, n.genes=ncol(M)) {
	stopifnot(ncol(M) > 1, n.genes <= ncol(M))
	if (!isSymmetric(M)) warning("Non-symmetric matrix: check your input, results are probably meaningless.")
	M <- M[1:n.genes, 1:n.genes]
	ee <- try(eigen(M), silent=TRUE)
	ee$vectors <- t (t(ee$vectors * sign(ee$vectors[2,]))) # by convention, the second eigenvector is always positive
	
	if (class(ee) == "try-error") {
		return(list(cor=NA, angle=NA, size=NA, eccentricity=NA))
	}
	ans <- list()
	ans$cor <- cov2cor(M)[upper.tri(M)]
	ans$angle <- modulopi(acos(ee$vectors[1,]))
	ans$size <- sum(diag(M))
	ans$eccentricity <- ee$values[2:ncol(M)]/ee$values[1]
	ans
}

 # Analytic methods
matrix2.noangle <- function(r, S, e) {
	cc <- if(r==0) 0 else r*S*sqrt(1-(1-e)^2/(1+e)^2)/(2*sqrt(1-r^2))
	v2 <- if (r==0) 
				e*S/(e+1)
			else
				(S+sqrt(S^2-4*cc^2/r^2))/2
	v1 <- S-v2
	matrix(c(v1, cc, cc, v2), ncol=2)
}

 # Numeric methods
matrix2.targetdist <- function(pp, target) {
	if (pp["v1"] <= 0 || pp["v2"] <= 0) return(Inf)
	mf <- matrix.features(matrix(c(pp["v1"], pp["c"], pp["c"], pp["v2"]), ncol=2))
	ans <- sum(sapply(names(target), function(nn) mf[[nn]][1]-target[nn])^2)
	ans
}

matrix2.optim <- function(target) {
	stopifnot(all(names(target) %in% c("cor","angle","size","eccentricity")))
	stopifnot(length(target) == 3)

	st <- c(target["size"]/2, target["size"]/2, sign(target["angle"])*target["size"]/10)
	names(st) <- c("v1","v2","c")
	ans <- optimx::optimx(par=st, matrix2.targetdist, method="Nelder-Mead", target=target, itnmax=2000)
	if (ans$convcode != 0) 
	  warning("Convergence failed for ", paste0(names(target), "=", target, collapse="  "))
	matrix(c(ans$v1, ans$c, ans$c, ans$v2), ncol=2)
}

##
 # Common wrapper function
matrix2.from.features <- function(cor=NA, angle=NA, size=NA, eccentricity=NA, method=c("numeric", "analytic")[1]) {
	if (method == "numeric") {
		tt <- c(cor=cor, angle=modulopi(angle), size=size, eccentricity=eccentricity)
		tt <- tt[!is.na(tt)]
		return(matrix2.optim(tt))
	} else if (method == "analytic") {
		if (is.na(cor) || is.na(size) || is.na(eccentricity) ||!is.na(angle))
			stop("Analytic method not available for this combination.")
		return(matrix2.noangle(r=cor, S=size, e=eccentricity))
	} else {
		stop("Method ", method, " not available.")
	}
}

##
#Function to print parameter files
param.from.sel.features <- function(param.template, param.out="param.par", cor=NA, angle=NA, size=NA, eccentricity=NA) {
	pp <- read.param(param.template)
	S.mat <- matrix2.from.features(cor=cor, angle=angle, size=size, eccentricity=eccentricity)
	n.genes <- if ("GENET_NBPHEN" %in% names(pp)) pp$GENET_NBPHEN else pp$GENET_NBLOC
	pS <- param.S.matrix(S.mat, n.genes=n.genes)
	pp[names(pS)] <- pS
	write.param(pp, parfile=param.out)
}


#########Function to draw ellipses axes###############
ellipse.axes <- function(x, scale=c(1,1), centre=c(0,0), level=0.95, t=sqrt(qchisq(level, 2)), which = c(1, 2)) {
  # Not sure that the scale factor works
  ee <- eigen(x[which, which])
  vvv <- t(t(ee$vectors)*sqrt(abs(ee$values))*t*scale)
  rbind(rbind(centre+vvv[,1], centre-vvv[,1]), c(NA,NA), rbind(centre+vvv[,2], centre-vvv[,2]))
}


######### Simplify fractions (to draw the angle axis in fractions of pi ############
######### From MASS package, copyright (C) 1994-2005 W. N. Venables and B. D. Ripley, GPL licence

fraction <- function(x, cycles = 10, max.denominator = 2000)
{
  a0 <- rep(0, length(x))
  A <- matrix(b0 <- rep(1, length(x)))
  fin <- is.finite(x)
  B <- matrix(floor(x))
  r <- as.vector(x) - drop(B)
  len <- 0
  while(any(which <- fin & (r > 1/max.denominator)) &&
	(len <- len + 1) <= cycles) {
    a <- a0
    b <- b0
    a[which] <- 1
    r[which] <- 1/r[which]
    b[which] <- floor(r[which])
    r[which] <- r[which] - b[which]
    A <- cbind(A, a)
    B <- cbind(B, b)
  }
  pq1 <- cbind(b0, a0)
  pq <- cbind(B[, 1], b0)
  len <- 1
  while((len <- len + 1) <= ncol(B)) {
    pq0 <- pq1
    pq1 <- pq
    pq <- B[, len] * pq1 + A[, len] * pq0
  }
  pq[!fin, 1] <- x[!fin]
  pq
}
