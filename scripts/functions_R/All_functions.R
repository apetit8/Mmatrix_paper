#FUNCTIONS

#Packages
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
library(Rcpp)
library(tidyr)
library(rlist)
library(tidyverse)
library(RColorBrewer)
library(data.table)

##TOOLS#######################################################################################

replicate.mean <- function(tt) {
  stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
  arr <- do.call(abind, c(tt, list(along=3)))
  rowMeans(arr, dims=2)
}

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

#Matrix extraction from output files
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

extract.fitness <- function(tt, what="Mfit", gen=tt[nrow(tt), "Gen"]) {
  ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
  if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
  return(ex)
}

extract.correlation <- function(tt, gen=tt[nrow(tt), "Gen"]) {
  M.mat <- extract.M.matrix(tt, gen=gen)
  ex <- (M.mat[1,2])/(sqrt(M.mat[1,1]) * sqrt(M.mat[2,2]) )
  return(ex)
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
  param.grid[,2] <- as.character(param.grid[,2] ) #the number = number of the argument : nextpar
  names.grid <- do.call(expand.grid,  list(lapply(param.list, function(pp) if (is.null(names(pp))) as.character(pp) else  as.character(pp))))
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
    # parfile <- paste0(param.dir, "/","simu-", rownames(param.grid)[i], param.ext)
    parfile <- paste0(param.dir, "/","simu-", str_split(rownames(param.grid)[i],"../", n=2, simplify = TRUE)[2], param.ext)
    write.param(parfile=parfile, parlist=mypar)
    for (rr in seq_len(reps)) {
      launch <- c(launch, paste(simevolv, "-p", paste0(rownames(param.grid)[i], param.ext), "-o", paste0(rownames(param.grid)[i], sep[3], formatC(rr, width = 1+floor(log10(reps)), format = "d", flag = "0"), simu.ext), sep=" "))
    }
  }
  # browser()
  if (!is.na(launchfile)) {
    dir <- as.character(dirname(template))
    cat(launch, file=print(sprintf("%s/%s", dir , launchfile)), sep="\n")
  }
  launch
}

#Function to print parameter files of the wanted S matrix
param.from.sel.features <- function(param.template, param.out="param.par", cor=NA, angle=NA, size=NA, eccentricity=NA, reps=40) {
  pp <- read.param(param.template)
  S.mat <- matrix2.from.features(cor=cor, angle=angle, size=size, eccentricity=eccentricity)
  n.genes <- if ("GENET_NBPHEN" %in% names(pp)) pp$GENET_NBPHEN else pp$GENET_NBLOC
  pS <- param.S.matrix(S.mat, n.genes=n.genes)
  pp[names(pS)] <- pS
  write.param(pp, parfile=param.out)
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

extract.S.matrix.fig3 <- function(parfile) {
  .cor2cov <- function(cc, sd) t(cc*sd)*sd
  rp <- read.param(parfile)
  fs <- rp$FITNESS_STRENGTH
  np <- 6 #manually = number of genes
    vv  <- 1/(2*fs)
    rr <- rp$FITNESS_CORRELATION
    ansr <- diag(np)
    ansr[upper.tri(ansr)] <- rr
    ansr <- as.matrix(Matrix::forceSymmetric(ansr))
    ans <- .cor2cov(ansr, sqrt(vv))
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
    #diag(dg) <- Inf
    diag(dg) <- Inf
    dg[1] <- 10 #To keep some selection on C and D genes, value can be changed
    dg[6] <- 10 #To keep some selection on C and D genes
    S.mat <- rbind(cbind(S.mat, rec), cbind(t(rec), dg))
  }
  list(
    FITNESS_STRENGTH = 1/(2*diag(S.mat)),
    FITNESS_CORRELATION = cov2cor(S.mat)[upper.tri(S.mat)])
}

#MATRIX FEATURES############################################################################
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

#Functions to find a matrix that gives an ellipse of given parameters (Numeric methods)
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

#MUTSD CALIBRATION###############################################################################
library(parallel)
sim.path <- "../../../simevolv/bin/Release/Simul_Prog"
can.tests <- 10

M.from.mutsd <- function(param.multilin, mutsd, gen=2000) {
  pp <- read.param(param.multilin)
  stopifnot(pp$TYPE_ARCHI == "multilinear")
  ff <- tempfile()
  pff <- paste0(ff, ".par")
  off <- paste0(ff, ".out")
  pp$SIMUL_GENER     <- gen
  pp$SIMUL_OUTPUT    <- gen
  pp$GENET_MUTSD     <- mutsd
  pp$OUT_CANAL_TESTS <- can.tests
  pp$OUT_CANAL_MUTSD <- mutsd
  write.param(pp, parfile=pff)
  system(paste(sim.path, "-p", pff, "-o", off, sep=" "))
  res <- read.table(off, header=TRUE)
  extract.M.matrix(res)
}

M.from.wagner <- function(param.wagner, gen=2000) {
  pp <- read.param(param.wagner)
  stopifnot(pp$TYPE_ARCHI %in% c("m2","wagner"))
  ff <- tempfile()
  pff <- paste0(ff, ".par")
  off <- paste0(ff, ".out")
  pp$SIMUL_GENER     <- gen
  pp$SIMUL_OUTPUT    <- gen
  pp$OUT_CANAL_TESTS <- can.tests
  write.param(pp, parfile=pff)
  system(paste(sim.path, "-p", pff, "-o", off, sep=" "))
  res <- read.table(off, header=TRUE)
  extract.M.matrix(res)	
}

M.from.multi <- function(param.multilin, gen=2000) {
  pp <- read.param(param.multilin)
  stopifnot(pp$TYPE_ARCHI %in% c("multilinear"))
  ff <- tempfile()
  pff <- paste0(ff, ".par")
  off <- paste0(ff, ".out")
  pp$SIMUL_GENER     <- gen
  pp$SIMUL_OUTPUT    <- gen
  pp$OUT_CANAL_TESTS <- can.tests
  write.param(pp, parfile=pff)
  system(paste(sim.path, "-p", pff, "-o", off, sep=" "))
  res <- read.table(off, header=TRUE)
  extract.M.matrix(res)	
}

mutsd.predict <- function(param.multilin, Vm, reps=10, gen=2000, mlim=c(0.0001, 1), traits=NA, mc.cores=detectCores()-1) {
  if (is.matrix(Vm)) Vm     <- diag(Vm)
  if (is.na(traits[1])) traits <- seq_along(Vm)
  mm <- exp(seq(log(mlim[1]), log(mlim[2]), length.out=reps))
  vm <- unlist(mclapply(mm, function(msd) {
    sum(diag(M.from.mutsd(param.multilin, msd, gen=gen))[traits])
  }, mc.cores=mc.cores))
  log.mutsd <- log(mm)
  log.Vm    <- log(vm)
  bb <- which.max(vm[vm<sum(Vm)])
  if (bb==1) bb <- 2
  if (bb==reps) bb <- reps-1
  log.mutsd <- log.mutsd[(bb-1):(bb+1)]
  log.Vm    <- log.Vm[(bb-1):(bb+1)]
  ll <- lm(log.mutsd ~ log.Vm)
  exp(predict(ll, list(log.Vm=log(sum(Vm)))))
}

calibrate.multilin <- function(param.multilin, param.wagner, gen=2000, traits=1:2, ...) {
  M.wagner <- M.from.wagner(param.wagner, gen=gen)
  Vm <- M.wagner[traits,traits]
  mst <- mutsd.predict(param.multilin, Vm, gen=gen, traits=traits, ...)
  pp <- read.param(param.multilin)
  pp$GENET_MUTSD     <- mst
  pp$OUT_CANAL_MUTSD <- mst
  write.param(pp, parfile=param.multilin)
}

#TRIGONOMETRY##########################################################################

# returns the angle of an ellipse (between -pi/2 and pi/2)
modulopi <- function(angle) {
  ans <- angle %% pi
  ifelse(ans > pi/2, ans-pi, ans)
}

#Diff between angle, allows to pick the modulo
modulo.all <- function(angle, modulo=pi) {
  angle <- angle %% modulo
  ifelse(angle > modulo/2, angle - modulo, angle)
}

#Calculate an angular mean between pi and -pi
mean.angle.2pi <- function(data) {
  atan2(mean(sin(data)), mean(cos(data)))
}

#Calculate an angular mean between pi/2 and -pi/2
mean.angle.pi <- function(data) {
  modulo.all(mean.angle.2pi(2*(data %% pi))/2, pi)
}

#Put mean M direction close to the S defined between -pi/2 and pi/2
mean.angle.pi.byS <- function(data, ang_S) {
  df <- modulo.all(mean.angle.2pi(2*(data %% pi))/2, pi)
  if ((df - mean(ang_S, 7)/2) < -pi/2){ df = df + pi}
  if ((df - mean(ang_S, 7)/2) > pi/2){ df = df - pi}
  return(df)
} 

#GRAPHIC FUNCTIONS###############################################################################
makeTransparent<-function(someColor, alpha=30)
{
  # From https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
  # Author: Nick Sabbe
  # Licence : CC-attribution-SA from the conditions of the website
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#Function to draw ellipses axes
ellipse.axes <- function(x, scale=c(1,1), centre=c(0,0), level=0.95, t=sqrt(qchisq(level, 2)), which = c(1, 2)) {
  # Not sure that the scale factor works
  ee <- eigen(x[which, which])
  vvv <- t(t(ee$vectors)*sqrt(abs(ee$values))*t*scale)
  rbind(rbind(centre+vvv[,1], centre-vvv[,1]), c(NA,NA), rbind(centre+vvv[,2], centre-vvv[,2]))
}

draw.ellipses <- function(G.mat=NULL, M.mat=NULL, S.mat=NULL, G.centre=c(0,0), M.centre=G.centre, S.centre=c(0,0),
                          M.factor=1, S.factor=0.00005, G.factor=1,
                          mat.col=c(G="dodgerblue4", M="darkolivegreen4", S="orange"), 
                          xlab="Trait 1", ylab="Trait 2", xlim=NULL, ylim=NULL, add=FALSE, asp=1, ...) { 
  if (!add) {
    if (is.null(xlim))
      xlim <- 2*c(-1,1)*2*sqrt(G.mat[1,1])
    if (is.null(ylim))
      ylim <- 2*c(-1,1)*2*sqrt(G.mat[2,2])
    plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, asp=asp, ...)
  }
  if (!is.null(G.mat) && all(diag(G.mat) > 1e-16)) {
    lines(ellipse::ellipse(G.factor*G.mat, centre=G.centre), col=mat.col["G"])
    lines(ellipse.axes(G.factor*G.mat, centre=G.centre), col=mat.col["G"])
  }
  if (!is.null(M.mat)) {
    lines(ellipse::ellipse(M.factor*M.mat, centre=M.centre), col=mat.col["M"])
    lines(ellipse.axes(M.factor*M.mat, centre=M.centre), col=mat.col["M"])
  }
  if (!is.null(S.mat)) {
    lines(ellipse::ellipse(S.factor*S.mat, centre=S.centre), col=mat.col["S"])
    lines(ellipse.axes(S.factor*S.mat, centre=S.centre), col=mat.col["S"])
  }
}

#Draw ellipses from multiple files into one plot
oneplot.allellipse <- function(data.dirs, xlim=NULL, ylim=NULL, xlab="Trait A", ylab="Trait B", Sell=TRUE, Mell=TRUE,
                               S.factor=0.00005, M.factor=1,  G.factor=1, xcoord=c(2000, 10000),ycoord=c(-1, 1.6),
                               asp=1, all.gen=FALSE, all.reps=FALSE, legend=TRUE, another_plot=FALSE, ...) {
  
  if (another_plot==TRUE){
  par(fig = c(grconvertX(xcoord, from="user", to="ndc"),
              grconvertY(ycoord, from="user", to="ndc")),
      mar = c(4,6,1,1),
      new = TRUE) #Ne fonctionne pas car : l'autre plot est tracé dans la fonction ... ...
    }
  
  plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, asp=asp,  ...)

  for (data.dir in data.dirs) {
    
    data.files = list.files(path=data.dir, pattern="simul.*.txt$", full.names=TRUE)
    param.file = list.files(path=data.dir, pattern="*.par",        full.names=TRUE)
    #browser()
    rp <- read.param(param.file)
    if (rp$TYPE_ARCHI %in% c("additive", "multilinear")) {
      mat.col=c(M="darkblue", S="orange") }
    if (rp$TYPE_ARCHI %in% c("m2")) {
      mat.col=c(M="yellowgreen", S="orange") }
    
    stopifnot(length(data.files) >0, length(param.file) == 1)
    
    simuls = lapply(data.files, function(x) read.table(x, header=TRUE)) #import list
    simuls.mean = replicate.mean(simuls)
    
    mygens <-rev(if (all.gen) simuls.mean[,"Gen"] else simuls.mean[nrow(simuls.mean),"Gen"])
    
    for (gen in mygens) {
      if (all.reps) {
        first <- TRUE
        for (simul in simuls) {
          phen.mean <- extract.P.mean(simul, gen=gen)
          
          draw.ellipses(
            M.mat=extract.M.matrix(simul, gen=gen),
            S.mat=NULL, G.centre=phen.mean, M.factor=M.factor, S.factor=S.factor,
            mat.col=makeTransparent(mat.col), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
            add=TRUE, ...)
          first <- FALSE          
        }
      }
      # Mean over replicates
      phen.mean <- extract.P.mean(simuls.mean, gen=gen)
      draw.ellipses(
        M.mat=if(Mell==TRUE) extract.M.matrix(simuls.mean, gen=gen) else NULL,
        S.mat=if(gen==mygens[1] && Sell==TRUE) extract.S.matrix(param.file) else NULL,
        G.centre=phen.mean, S.centre=extract.theta(param.file),
        M.factor=M.factor, S.factor=S.factor,
        mat.col=mat.col, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
        add=TRUE, ...)
      if (legend==TRUE){legend("topleft", lty=1, box.lty=0, bg="transparent", col=c("yellowgreen", "darkblue", "orange"),
             legend=c(paste0("M GRN", if (M.factor != 1) paste0(" (x ", M.factor, ")")), paste0("M multilinear", if (M.factor != 1) paste0(" (x ", M.factor, ")")),paste0("S", if (S.factor != 1) paste0(" (x ", S.factor, ")"))))}
    }
  }
}

oneplot.netw <- function(data.dirs, xlim=NULL, ylim=NULL, 
                               M.factor=1,  G.factor=1, S.factor=1,
                               asp=1, all.gen=FALSE, all.reps=TRUE, mat.col=c(M="yellowgreen", S="orange"), angle=c(1,0), ...) {
  for (data.d in data.dirs) {
    for (i in angle) {
    # add <- FALSE
    data.dir  <- list.files(path=data.d, pattern=sprintf("angle%s.parangle", i), full.names=TRUE)
    print(data.dir)
    data.files <- list.files(path=data.dir, pattern="simul.*.txt$", full.names=TRUE)
    data.files <- data.files[sapply(data.files, file.size) > 10000]
    param.file <- list.files(path=data.dir, pattern="*.par",        full.names=TRUE)
    # browser()
    stopifnot(length(data.files) >0, length(param.file) == 1)
    
    simuls = lapply(data.files, function(x) read.table(x, header=TRUE)) #import list
    simuls.mean = replicate.mean(simuls)
    
    mygens <-rev(if (all.gen) simuls.mean[,"Gen"] else simuls.mean[nrow(simuls.mean),"Gen"])
    plot(NULL,  xlim=xlim, ylim=ylim, asp=asp, ...)
    
    for (gen in mygens) {
      # Mean over replicates
      phen.mean <- extract.P.mean(simuls.mean, gen=gen)
      draw.ellipses(
        M.mat=NULL,
        S.mat=if(gen==mygens[1]) extract.S.matrix(param.file) else NULL,
        G.centre=extract.theta(param.file), S.centre=extract.theta(param.file),
        M.factor=M.factor, S.factor=S.factor,
        mat.col=mat.col, xlim=xlim, ylim=ylim,
        add=TRUE, ...)
      if (all.reps) {
        for (simul in simuls) {
          phen.mean <- extract.P.mean(simul, gen=gen)
          # browser()
          draw.ellipses(
            M.mat=extract.M.matrix(simul, gen=gen),
            S.mat=NULL, G.centre=extract.theta(param.file), M.factor=M.factor, S.factor=S.factor,
            mat.col=makeTransparent(mat.col, alpha=40), xlim=xlim, ylim=ylim,
            add=TRUE, ...)
          # add <- TRUE         
        }
      }
    }
   }
  }
}


##DATAFRAME#########################################################
diff.ms.df <- function(sims.dir, modulo=pi){
  simul.df <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
    files <- files[sapply(files, file.size) > 500000] #10000
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
    files <- files[sapply(files, file.size) > 10000]
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
# 
# df.raster.map <- function(sims.dir, modulo=pi, all.gen=FALSE){
#   simul.df <- data.frame()
#   filedata <- data.frame()
#   #collect all the data
#     files     <- list.files(path = sims.dir, full.names=TRUE, pattern = "\\.txt$")
#     #change here how to get to the right param file !
#     files <- files[sapply(files, file.size) > 150000]
#     for (i in files) { 
#         #
#         parname1<- str_split(i, ".txt$", n=2, simplify = TRUE)
#         parname2<- str_split(parname1, "-R", n=2, simplify = TRUE)
#         parname3<- str_split(parname2, "NEXTPAR-", n=2, simplify = TRUE)
#         parname4 <- str_split(parname2, "/FITNESS_OPTIMUM-", n=2, simplify = TRUE)
#         param.file <- sprintf("%s/%s", parname4[1,1], parname3[1,2]) #get the name of the second paramfile, with the right S
#         #
#         rp <- read.param(param.file)
#         S.ref2 <- matrix.features(extract.S.matrix.fig3(param.file), n.genes=2)[["angle"]][1]
#         
#         mytopos <- lapply(i, function(ff) {
#           tt <- read.table(ff, header=TRUE) 
#           mygens <-rev(if (all.gen==TRUE) tt[,"Gen"] else tt[nrow(tt),"Gen"])
#           for (gen in mygens) {
#             phen.mean <- extract.P.mean(tt, gen=gen)
#             M.mat <- extract.M.matrix(tt, gen=gen)
#             M.feat2 <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M angle
#             mean.ma.df <- M.feat2
#             data.gen <- c(ff, gen, S.ref2, M.feat2, phen.mean[1],phen.mean[2],mean.ma.df )
#             filedata <- rbind(filedata, data.gen)
#             #create a list of data
#           }
#           return(as.list(filedata))
#     })
#     newt <- as.data.frame(rbindlist(mytopos, use.names=FALSE))
#     newt[,7] <- mean.angle.pi(as.numeric(newt[,7]))
#     setnames(newt, 1:7, c("data.dir", "Gen","ang_S","ang_M","P_mean_A","P_mean_B","mean_ang_m"))
#     # newt <- ldply(mytopos) # create a data.frame for mytopos
#     simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
#     }
#   simul.df[,2:7] <- lapply( simul.df[,2:7], as.numeric)
#   sq_dist <- 1 - (((modulo.all(simul.df$ang_M-simul.df$ang_S, modulo = modulo))^2)/((pi^2)/12))
#   #sq_dist <- modulo.all(simul.df$ang_M-simul.df$ang_S, modulo = modulo)
#   simul.df <- as.data.frame(cbind(simul.df, sq_dist)) #df of 8 columns
#   return(simul.df)
# }

df.raster.map <- function(sims.dir, modulo=pi, all.gen=FALSE){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  files     <- list.files(path = sims.dir, full.names=TRUE, pattern = "\\.txt$")
  #change here how to get to the right param file !
  files <- files[sapply(files, file.size) > 150000]
  for (i in files) { 
    #
    parname1<- str_split(i, ".txt$", n=2, simplify = TRUE)
    parname2<- str_split(parname1, "-R", n=2, simplify = TRUE)
    parname3<- str_split(parname2, "NEXTPAR-", n=2, simplify = TRUE)
    parname4 <- str_split(parname2, "/FITNESS_OPTIMUM-", n=2, simplify = TRUE)
    param.file <- sprintf("%s/%s", parname4[1,1], parname3[1,2]) #get the name of the second paramfile, with the right S
    #
    rp <- read.param(param.file)
    S.ref2 <- matrix.features(extract.S.matrix.fig3(param.file), n.genes=2)[["angle"]][1]
    
    mytopos <- lapply(i, function(ff) {
      tt <- read.table(ff, header=TRUE) 
      mygens <-rev(if (all.gen==TRUE) tt[,"Gen"] else tt[nrow(tt),"Gen"])
      for (gen in mygens) {
        phen.mean <- extract.P.mean(tt, gen=gen)
        M.mat <- extract.M.matrix(tt, gen=gen)
        M.feat2 <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M angle
        fitness <- extract.fitness(tt, "MFit", gen=gen)
        data.gen <- c(ff, gen, S.ref2, M.feat2, phen.mean[1],phen.mean[2], fitness)
        filedata <- rbind(filedata, data.gen)
        #create a list of data
      }
      return(as.list(filedata))
    })
    newt <- as.data.frame(rbindlist(mytopos, use.names=FALSE))
    #newt[,7] <- mean.angle.pi(as.numeric(newt[,7]))
    setnames(newt, 1:7, c("data.dir", "Gen","ang_S","ang_M","P_mean_A","P_mean_B","Fitness"))
    # newt <- ldply(mytopos) # create a data.frame for mytopos
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  simul.df[,2:7] <- lapply( simul.df[,2:7], as.numeric)
  sq_dist <- 1 - (((modulo.all(simul.df$ang_M-simul.df$ang_S, modulo = modulo))^2)/((pi^2)/12))
  # sq_dist <- modulo.all(simul.df$ang_M-simul.df$ang_S, modulo = modulo)
  simul.df <- as.data.frame(cbind(simul.df, sq_dist)) #df of 8 columns
  return(simul.df)
}

df.fig1 <- function(sims.dir, modulo=pi, all.gen=FALSE){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    print(i)
    files     <- list.files(path = i, full.names=TRUE, pattern = "\\.txt$")
    files <- files[sapply(files, file.size) > 1000]
    param.file.all = list.files(path=i, pattern="\\.par$", full.names=TRUE)
    param.file <- as.character(param.file.all[1])
    stopifnot(length(files) >0)
    
    #read param.par to associate each simu with S properties
    rp <- read.param(param.file)
    mut.rate <- 2*rp$GENET_MUTRATES
    if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    S.ref2 <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["angle"]][1]
    
    mytopos <- lapply(files, function(ff) {
      #print(ff)
      tt <- read.table(ff, header=TRUE) 
      mygens <-rev(if (all.gen==TRUE) tt[,"Gen"] else tt[nrow(tt),"Gen"])
      for (gen in mygens) {
        phen.mean <- extract.P.mean(tt, gen=gen)
        M.mat <- extract.M.matrix(tt, gen=gen)
        M.feat2 <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M angle
        mean.ma.df <- M.feat2
        data.gen <- c(ff, gen, S.ref2, M.feat2, phen.mean[1],phen.mean[2],mean.ma.df )
        filedata <- rbind(filedata, data.gen)
        #create a list of data
      }
      return(as.list(filedata))
    })
    newt <- as.data.frame(rbindlist(mytopos, use.names=FALSE))
    newt[,7] <- mean.angle.pi(as.numeric(newt[,7]))
    setnames(newt, 1:7, c("data.dir", "Gen","ang_S","ang_M","P_mean_A","P_mean_B","mean_ang_m"))
    # newt <- ldply(mytopos) # create a data.frame for mytopos
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  simul.df[,2:7] <- lapply( simul.df[,2:7], as.numeric)
  sq_dist <- 1 - (((modulo.all(simul.df$ang_M-simul.df$ang_S, modulo = modulo))^2)/((pi^2)/12))
  simul.df <- as.data.frame(cbind(simul.df, sq_dist)) #df of 8 columns
  return(simul.df)
}


df.topo.wide.m<- function(sims.dir, w_of_4=FALSE, w_of_6=FALSE, network=FALSE, file_size=100000){
  simul.df <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    print(i)
    files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
    files <- files[sapply(files, file.size) > file_size]
    param.file = list.files(path=i, pattern="param.par", full.names=TRUE)
    # stopifnot(length(files) >0, length(param.file) == 1)
    #read param.par to associate each simu with the S "what"
    rp <- read.param(param.file)
    mut.rate <- 2*rp$GENET_MUTRATES
    if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    S.ref2 <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["angle"]][1] #S matrix
    
    mytopos <- lapply(files, function(ff) {
      print(ff)
      tt <- read.table(ff, header=TRUE)
      #M data
      gen <-tt[nrow(tt),"Gen"]
      M.mat <- extract.M.matrix(tt, gen=gen)
      M.ang_ppi <- matrix.features(M.mat,n.genes=2)[["angle"]][1]+pi #M matrix
      M.ang_mpi <- matrix.features(M.mat,n.genes=2)[["angle"]][1]-pi #M matrix
      M.feat1 <- matrix.features(M.mat,n.genes=2)[["eccentricity"]][1] #M matrix
      M.ang <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
      M.feat3 <- matrix.features(M.mat,n.genes=2)[["size"]][1] #M matrix
      fitness <- extract.fitness(tt, "MFit", gen=gen)
      # opt1 <- extract.fitness(tt, "FitOpt1", gen=gen)[1]
      corr <- extract.correlation(tt, gen=gen)
      opt2 <- extract.fitness(tt, "FitOpt2", gen=gen)[1]
      if (network==TRUE){
        Wmat <- extract.W.matrix(tt) #W data
        new <- c(ff, M.ang_ppi, S.ref2, M.ang_mpi, M.feat1, M.ang, M.feat3, fitness, corr, Wmat )
      }
      if (network!=TRUE){
        new <- c(ff, M.ang_ppi, S.ref2, M.ang_mpi, M.feat1, M.ang, M.feat3, fitness, corr )
      }
      return(new)
      #create a list of data
    })
    newt <- ldply(mytopos) # create a data.frame for i
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  
  if (network!=TRUE){
    simul.df <- setnames(simul.df[,1:9], c("data.dir","ang_M_ppi","ang_S","ang_M_mpi","ecc_M","ang_M","siz_M","fitness","corr"))
    simul.df[,2:9] <- lapply( simul.df[,2:9], as.numeric)
  }
  if (w_of_6==TRUE & network==TRUE) {
    simul.df <- setnames(simul.df[,1:45], c("data.dir","ang_M_ppi","ang_S","ang_M_mpi","ecc_M","ang_M","siz_M","fitness","corr",
                               "A_A","B_A","C_A","D_A","E_A","F_A","A_B","B_B","C_B","D_B","E_B","F_B","A_C","B_C","C_C","D_C","E_C","F_C",
                               "A_D","B_D","C_D","D_D","E_D","F_D","A_E","B_E","C_E","D_E","E_E","F_E","A_F","B_F","C_F","D_F","E_F","F_F"))
    simul.df[,2:45] <- lapply( simul.df[,2:45], as.numeric)
  }
  if (w_of_4==TRUE & network==TRUE) {
    #Naming network cells (for 4 genes network)
    names(simul.df)[names(simul.df) == "V10"] <- "B_A"
    names(simul.df)[names(simul.df) == "V11"] <- "C_A"
    names(simul.df)[names(simul.df) == "V12"] <- "D_A"
    names(simul.df)[names(simul.df) == "V13"] <- "A_B"
    names(simul.df)[names(simul.df) == "V15"] <- "C_B"
    names(simul.df)[names(simul.df) == "V16"] <- "D_B"
    names(simul.df)[names(simul.df) == "V17"] <- "A_C"
    names(simul.df)[names(simul.df) == "V18"] <- "B_C"
    names(simul.df)[names(simul.df) == "V20"] <- "D_C"
    names(simul.df)[names(simul.df) == "V21"] <- "A_D"
    names(simul.df)[names(simul.df) == "V22"] <- "B_D"
    names(simul.df)[names(simul.df) == "V23"] <- "C_D"
    simul.df[,9:24] <- lapply( simul.df[,9:24], as.numeric)
  }
  return(simul.df)
}

#Function for optimum map alignment data
df.opt.map <- function(sims.dir, modulo=pi, bymean=FALSE){
  df.all <- data.frame(NULL)
  for (i in sims.dir) {
    df <- df.fig1(i, modulo=modulo)
    pop <- str_split(i, "../../simul/fig_3/r", n=2, simplify = TRUE)
    df[,9] <- sprintf("%s", pop[,2])
    opti <- sub("-R.*", "", df$data.dir)
    df[,10] <- sprintf("%s", opti)
    df.all <- rbind(df.all, df)
  }
  return(df.all)
  #This part add a column : Xi_a calculated from the mean M angle. 
  if (bymean==TRUE){
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
}


#Function for optimum map alignment data
df.opt.map2 <- function(sims.dir, modulo=pi, bymean=FALSE, gen=FALSE){
  df.all <- data.frame(NULL)
    df <- df.raster.map(sims.dir, modulo=modulo, all.gen=gen)
    pop <- str_split(sims.dir, "../../simul/fig_3_bis/", n=2, simplify = TRUE)
    df$pop <- sprintf("%s", pop[,2])
    opti <- sub("-R.*", "", df$data.dir)
    df$opt <- sprintf("%s", opti)
    df.all <- rbind(df.all, df)
  return(df.all)
  #This part add a column : Xi_a calculated from the mean M angle. 
  if (bymean==TRUE){
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
}

df.data <- function(sims.dirs, modulo=pi, pattern="../../simul/", variable="population", network=FALSE, w_of_6=FALSE, file_size=100000){
  df.datas <- data.frame(NULL)
  for (i in sims.dirs) {
    sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
    df <- df.topo.wide.m(sims.dir, network=network, w_of_6=w_of_6, file_size=file_size)
    # df[,9] <- modulo.all(df$ang_M-df$ang_S, modulo = modulo)
    # sq_rho <- round(mse.ms(df[,9]), 4)
    pop <- str_split(i, pattern, n=2, simplify = TRUE)
    df$pop <- sprintf("%s", pop[,2])
    df.datas <- rbind(df.datas, df)
  }
  # names(df.datas)[names(df.datas) == "V9"] <- "ang_diff"
  # names(df.datas)[names(df.datas) == "V10"] <- variable
  return(df.datas)
}

table.index <- function(sims.dirs, modulo=pi, pattern="../../simul/", ref=df.all, bymean=FALSE, asgrob=TRUE){
  ##Alignment index calcul and Tab
  df.sq_rho <- data.frame(ncol(4))
  j <- 1
  modulo <- pi
  for (i in sims.dirs) {
    sims.dir <- list.dirs(i, recursive = FALSE)
    if (bymean==TRUE){
      diff <- diff.ms.by.mean(sims.dir=sims.dir, modulo = modulo)}
    if (bymean==FALSE){
      diff <- diff.ms.df(sims.dir=sims.dir, modulo = modulo)}
    sq_rho <-round( mse.ms(diff), 4)
    delta_diff <-mse.ms(diff)/((pi^2)/12) #distance if M angle is constant
    
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
  # if (asgrob==TRUE) g <- tableGrob(as.matrix(df.sq_rho),rows = NULL) else g <- df.sq_rho
  g <- if(asgrob==TRUE) tableGrob(as.matrix(df.sq_rho),rows = NULL) else df.sq_rho
  return(g)
}

mse.ms <- function(diff.ls){
  mse <- mean((diff.ls)^2)
  return(mse)
} 

#GGPLOT2 CUSTOM FUNCTIONS############################################################
# Done with the precious help of this webpage : https://eliocamp-github-io.translate.goog/codigo-r/2018/06/tu-propio-geom-smooth/?_x_tr_sl=auto&_x_tr_tl=en&_x_tr_hl=fr&_x_tr_pto=nui

#mean for "_smooth" ggplot2 functions
AngleSmooth <- function(formula, data, weights, n = 0.5) {
  ff <- data.frame()
  yi <- unique(round(data$x, 2))
  for (k in yi) {
    ww <- subset(data, round(x, 2) == k)
    f <- mean.angle.pi(data=ww$y)
    dd <- data.frame(matrix(ncol = 2, nrow = nrow(ww)))
    dd[1:nrow(ww),1] <- f
    dd[,2] <- ww$x
    setnames(dd, 1:2, c("y","x"))
    ff <- rbind(ff, dd)
    
  }
  model <- list(x = ff$x, pred = ff$y)
  class(model) <- "angle_smooth"
  return(model)
}

predictdf.angle_smooth <- function(model, xseq, se=FALSE, level) {
  data.frame(x = model$x, y = model$pred)
}

#mean angle for "summary" ggplot2 functions
Angle_Mean <- function(data) {
  data.frame(y = mean.angle.pi(data))
}

#standard deviation for error bars in stat_summary ggplot2 function 
st_dev_angle <- function(data){
  sd <- sqrt( mean(modulo.all((data - mean.angle.pi(data)))^2 ))
  data.frame( ymin=(mean.angle.pi(data))-sd,
              ymax=(mean.angle.pi(data))+sd
  )
}

#standard deviation value for stat_summary ggplot2 function 
st_dev_abs <- function(data){
  sd <- sqrt( mean(modulo.all((data - mean.angle.pi(data)))^2 ))
  data.frame( sda = sd
  )
}

