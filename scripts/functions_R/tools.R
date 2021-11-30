suppressPackageStartupMessages(library(abind))
library(plyr)
library(stringr)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(matrixStats)
library(ggplot2)
library(ggExtra)
library(grid)
library(ggstance)
library(gtable)
library(scales)

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

#function to write multiple parameter files for multiple changing parameters
generate.param.list <- function(template, param.list, reps=1, launchfile=NA, sep=c("-", "-", "-R"), vec.indx=rep(1, length(param.list)), 
                                simevolv="$SIMEVOLV", param.dir=dirname(template), simu.dir=dirname(template), param.ext=".par", simu.ext=".txt") {
  
  launch <- NULL
  
  templ <- read.param(template)
  param.grid <- do.call(expand.grid, param.list)
  names.grid <- do.call(expand.grid, lapply(param.list, function(pp) if (is.null(names(pp))) as.numeric(pp) else names(pp)))
  colnames(names.grid) <- paste0(colnames(names.grid), ifelse(vec.indx == 1, "", vec.indx))
  rownames(param.grid) <- do.call(paste, c(lapply(colnames(names.grid), function(nn) paste(nn, names.grid[,nn], sep=sep[1])), list(sep=sep[2])))
  
  for (i in seq_len(nrow(param.grid))) {
    mypar <- templ
    for (j in seq_along(names(param.list))) {
      if (!names(param.list)[j] %in% names(mypar))
        stop("Error: parameter ", names(param.list)[j], " is not in file ", template, ".")
      if (length(mypar[[names(param.list)[j]]]) < vec.indx[j])
        stop("Error: parameter ", names(param.list)[j], " does not have ", vec.indx[j], "elements.")
      mypar[[names(param.list)[j]]][vec.indx[j]] <- param.grid[i,j]
    }
    parfile <- paste0(param.dir, "/", rownames(param.grid)[i], param.ext)
    write.param(parfile=parfile, parlist=mypar)
    for (rr in seq_len(reps)) {
      launch <- c(launch, paste(simevolv, "-p", paste0(rownames(param.grid)[i], param.ext), "-o", paste0(rownames(param.grid)[i], sep[3], formatC(rr, width = 1+floor(log10(reps)), format = "d", flag = "0"), simu.ext), sep=" "))
    }
  }
  if (!is.na(launchfile)) {
    dir <- as.character(dirname(template))
    cat(launch, file=print(sprintf("%s/%s", dir , launchfile)), sep="\n")
  }
  launch
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

######### Angles features

modulopi <- function(angle) {
	# returns the angle of an ellipse (between -pi/2 and pi/2)
	ans <- angle %% pi
	ifelse(ans > pi/2, ans-pi, ans)
}

modulo <- function(diff_angle, modulo=pi) {
  diff_angle <- diff_angle %% modulo
  ifelse(diff_angle > modulo/2, diff_angle - modulo, diff_angle)
}

modulo.all <- function(angle, modulo=pi) {
  angle <- angle %% modulo
  ifelse(angle > modulo/2, angle - modulo, angle)
}

mean.angle.2pi <- function(data) {
  atan2(mean(sin(data)), mean(cos(data)))
}

mean.angle.pi <- function(data) {
  modulo.all(mean.angle.2pi(2*(data %% pi))/2, pi)
}

angle.diff <- function(ang1, ang2) {
  diff <- ((pi/2)/(2*pi)) * acos( cos((2*pi*(ang1-ang2))/(pi/2)))
  return(diff)
}

######### Matrix features
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


diff.ms.df <- function(sims.dir, modulo=pi){
  simul.df <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
    param.file = list.files(path=i, pattern="param.par", full.names=TRUE)
    stopifnot(length(files) >0, length(param.file) == 1)
    
    #read param.par to associate each simu with the S "what"
    rp <- read.param(param.file)
    mut.rate <- 2*rp$GENET_MUTRATES
    if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    S.ref <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["angle"]][1]
    
    mytopos <- lapply(files, function(ff) {
      #print(ff)
      tt <- read.table(ff, header=TRUE)
      gen <-tt[nrow(tt),"Gen"]
      P.mat <- extract.P.matrix(tt, gen=gen)
      M.mat <- extract.M.matrix(tt, gen=gen)
      M.feat <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
      new <- c(S.ref, M.feat) #create a list of data
    })
    newt <- ldply(mytopos) # create a data.frame for i
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  diff.df <- data.frame()
  simul.df[,1:2] <- lapply( simul.df[,1:2], as.numeric)
  diff.ls <-modulo.all(simul.df[,2]-simul.df[,1], modulo=modulo)
  return(diff.ls)
}


diff.ms.by.mean <- function(sims.dir, modulo=pi){
  simul.df <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
    param.file = list.files(path=i, pattern="param.par", full.names=TRUE)
    stopifnot(length(files) >0, length(param.file) == 1)
    
    #read param.par to associate each simu with the S "what"
    rp <- read.param(param.file)
    mut.rate <- 2*rp$GENET_MUTRATES
    if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    S.ref <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["angle"]][1]
    
    mytopos <- lapply(files, function(ff) {
      #print(ff)
      tt <- read.table(ff, header=TRUE)
      gen <-tt[nrow(tt),"Gen"]
      P.mat <- extract.P.matrix(tt, gen=gen)
      M.mat <- extract.M.matrix(tt, gen=gen)
      M.feat <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
      empt <- M.feat
      new <- c(S.ref, M.feat, empt) #create a list of data
    })
    newt <- ldply(mytopos) # create a data.frame for i
    newt[,3] <- mean.angle.pi(newt[,3])
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  diff.df <- data.frame()
  simul.df[,1:3] <- lapply( simul.df[,1:3], as.numeric)
  diff.ls <-modulo.all(simul.df[,3]-simul.df[,1], modulo=modulo)
  return(diff.ls)
}

#Function for optimum map alignment data
df.opt.map <- function(sims.dir, modulo=pi){
  df.all <- data.frame(NULL)
  for (i in sims.dir) {
    df <- df.raster.map(i, modulo=modulo)
    pop <- str_split(i, "../../simul/fig_3/raster_initcorr", n=2, simplify = TRUE)
    df[,9] <- sprintf("%s", pop[,2])
    opti <- sub("-R.*", "", df$data.dir)
    df[,10] <- sprintf("%s", opti)
    df.all <- rbind(df.all, df)
  }
  
  df.mean.Ea <- data.frame()  
  for (j in unique(df.all$V9)) {
    df.ang <- subset(df.all, V9 == j)
    for (k in unique(df.ang$V10)) {
      df.ang.opt <- subset(df.ang, V10 == k)
      df.ang.opt[,11] <- 1 - ((modulo.all(mean.angle.pi(df.ang.opt$ang_M)-df.ang.opt$ang_S, modulo = modulo))^2)/((pi^2)/12)
      df.mean.Ea <- rbind(df.mean.Ea, df.ang.opt)
    }
  }
  return(df.mean.Ea)
}


table.index <- function(sims.dirs, modulo=pi, pattern="../../simul/", ref=df.all, bymean=TRUE){
    ##Alignment index calcul and Tab
    df.sq_rho <- data.frame(ncol(4))
    j <- 1
    modulo <- pi
    for (i in sims.dirs) {
      sims.dir <- list.dirs(i, recursive = FALSE)
      if (bymean==TRUE){
      diff <- diff.ms.by.mean(sims.dir=sims.dir, modulo = modulo)}
      if (bymean!=TRUE){
        diff <- diff.ms.df(sims.dir=sims.dir, modulo = modulo)}
      sq_rho <-round( mse.ms(diff), 4)
      delta_diff <-mean(mse.ms(diff)/(((pi^2)/12))) #distance if M angle is constant
      
      #mean square distance
      sdall <- data.frame()
      sddf <- filter(ref, grepl(i,sprintf("%s/", data.dir)))
      yi <- unique(round(sddf$ang_S, 3))
      for (l in yi) {
        ww <- subset(sddf, round(ang_S, 3) == l)
        sdy <- as.data.frame(sqrt( modulo.all((ww$ang_M - mean.angle.pi(ww$ang_M)))^2))
        sdall <- rbind(sdall, sdy)
      }
      sd <- mean(sdall[,1]) #square deviation
      
      pop <- str_split(i, pattern, n=2, simplify = TRUE)
      df.sq_rho[j,1] <- pop[2]
      df.sq_rho[j,2] <- sq_rho
      df.sq_rho[j,3] <- round(sd, 4)
      df.sq_rho[j,4] <- round(1 - delta_diff, 4)
      j <- j +1
    }
    names(df.sq_rho) <- c("Population","sq_rho","standard deviation","Xi_alpha")
    g <- tableGrob(as.matrix(df.sq_rho),rows = NULL)
    return(g)
}


df.data <- function(sims.dirs, modulo=pi, pattern="../../simul/", variable="population"){
    df.datas <- data.frame(NULL)
    for (i in sims.dirs) {
      sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
      df <- df.topo.raw(sims.dir, network=FALSE)
      df[,9] <- modulo.all(df$ang_M-df$ang_S, modulo = modulo)
      sq_rho <- round(mse.ms(df[,9]), 4)
      pop <- str_split(i, pattern, n=2, simplify = TRUE)
      df[,10] <- sprintf("%s", pop[,2])
      df.datas <- rbind(df.datas, df)
    }
    names(df.datas)[names(df.datas) == "V9"] <- "ang_diff"
    names(df.datas)[names(df.datas) == "V10"] <- variable
    return(df.datas)
}

#Calculate the MSE between M and S angle
mse.ms <- function(diff.ls){
  mse <- mean((diff.ls)^2)
  return(mse)
}

## Convert x-ticks to fractional x-ticks with a symbol multiplier
fracAx <- function(p, symbol=pi, y=FALSE, x=TRUE, width=0.25) {
  require(MASS)                                                 # for fractions
  val <- tryCatch(eval(parse(text=symbol)), error=function(e) 1)
  info <- ggplot_build(p)
  xrange <- range(pi/2,-pi/2)                 # get the x-range of figure
  vec.breaks <- seq(floor(xrange[1]), ceiling(xrange[2]), by=width)
  fracs <- strsplit(attr(fractions(vec.breaks), "fracs"), "/")  # convert to fractions
  if(y==TRUE){
  labels <- sapply(fracs, function(i)
    if (length(i) > 1) { paste(i[1], "*", symbol, "/", i[2]) }
    else { paste(i, "*", symbol) })
  if(x==TRUE){
    p <- p + scale_x_continuous(breaks=vec.breaks*val, labels=parse(text=labels))
    }
    p <- p + scale_y_continuous(breaks=vec.breaks*val, labels=parse(text=labels))
    return(p)
  }
  if(y!=TRUE){
    labels <- sapply(fracs, function(i)
      if (length(i) > 1) { paste(i[1], "*", symbol, "/", i[2]) }
      else { paste(i, "*", symbol) })
    p + scale_x_continuous(breaks=vec.breaks*val, labels=parse(text=labels))
  }
  
}