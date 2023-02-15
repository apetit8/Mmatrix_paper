#FUNCTIONS

#Packages
suppressPackageStartupMessages(library(abind))
library(plyr)
library(stringr) #
# library(FactoMineR)
# library(factoextra) 
library(gridExtra)
library(matrixStats) #
library(ggplot2) #
library(ggExtra)
library(grid)
# library(ggstance)
library(Rcpp)
library(tidyr)
library(rlist)
library(dplyr)
# library(tidyverse)
library(RColorBrewer) #
library(data.table) #

##TOOLS#######################################################################################


##PARAMFILE###################################################################################

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
param.from.sel.features <- function(param.template, param.out="param.par", cor=NA,
        angle=NA, size=NA, eccentricity=NA, reps=40, sel.genes=c(3:(2+(n.genes-2)/2)), initW=NULL ) {
  pp <- read.param(param.template)
  if(!is.null(initW)){pp$INIT_ALLELES_FULL <- initW}
  S.mat <- matrix2.from.features(cor=cor, angle=angle, size=size, eccentricity=eccentricity)
  n.genes <- if ("GENET_NBPHEN" %in% names(pp)) pp$GENET_NBPHEN else pp$GENET_NBLOC
  pS <- param.S.matrix(S.mat, n.genes=n.genes, sel.genes=sel.genes )
  pp[names(pS)] <- pS
  write.param(pp, parfile=param.out)
}


corparam <- function(nloc=6, mutsd=0.036865857, eccentricity=0.99) {
  angles <- seq(-pi/2+1e-6, pi/2-1e-6, length.out=nloc+1)[-1]
  template.mat <- 2*mutsd*mutsd*matrix2.from.features(angle=0, size=1, eccentricity=eccentricity)
  Ms <- lapply(angles, function(alpha) {
    rotation <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), ncol=2)
    rotation %*% template.mat %*% t(rotation)
  })
  ans <- c(
    paste0("GENET_MUTCOR     ", paste0(sapply(Ms, function(M) cov2cor(M)[1,2]), collapse=" ")),
    paste0("GENET_MUTSD      ", paste0(sapply(Ms, function(M) paste(sqrt(diag(M)), sep=" ")), collapse=" "))
  )
  return(ans)
}

##EXTRACT DATA####################################################################

replicate.mean <- function(tt) {
  stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
  arr <- do.call(abind, c(tt, list(along=3)))
  rowMeans(arr, dims=2)
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

extract.nbphen <- function(parfile) {
  rp <- read.param(parfile)
  if (rp$TYPE_ARCHI %in% c("additive", "multilinear","fkl")) {
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


#PHENOTYPE############################################################################
library(Rcpp)
library(inline, quietly=TRUE)

sigma.M2 <- function(x, a) {
  1. / (1. + exp((-x/(a*(a-1)))+log(1/a-1)))
}

sigma.M2p <- function(x, lambda, mu) {
  1. / (1. + lambda * exp(-mu*x))
}

suppressMessages(library(compiler))
sigma.M2c <- cmpfun(sigma.M2p)

#Take properties of S and convert it in S parameters for Simevolv 
param.S.matrix <- function(S.mat, n.genes=ncol(S.mat), sel.genes=c(3:(2+(n.genes-2)/2)) ) {
  stopifnot(is.matrix(S.mat), ncol(S.mat) == nrow(S.mat))
  stopifnot(n.genes > 0)
  if (n.genes < ncol(S.mat))
    S.mat <- S.mat[1:n.genes, 1:n.genes]
  if (n.genes > ncol(S.mat)) {
    rec <- matrix(0, ncol=n.genes-ncol(S.mat), nrow=nrow(S.mat))
    dg <- matrix(0, ncol=n.genes-ncol(S.mat), nrow=n.genes-nrow(S.mat))
    #diag(dg) <- Inf
    diag(dg) <- Inf
    S.mat <- rbind(cbind(S.mat, rec), cbind(t(rec), dg))
  }
  fit.strength <- 1/(2*diag(S.mat))
  fit.strength[sel.genes] <- 0.05
  list(
    FITNESS_STRENGTH = fit.strength,
    FITNESS_CORRELATION = cov2cor(S.mat)[upper.tri(S.mat)])
}

#R version of simevolv way to compute phenotype
internal_loop_R <- function(W, S0, a, steps, measure) {
  lambda <- (1-a)/a
  mu <- 1/(a*(1-a))
  sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
  sto[,1] <- S0
  for (i in 1:steps) {
    S0 <- sigma.M2c((W %*% S0), lambda=lambda, mu=mu) 			
    sto[,i+1] <- S0
  }
  list(full=sto, sumx=rowSums(sto[,(steps-measure+1):steps]), sumx2=rowSums(sto[,(steps-measure+1):steps]^2))
}

#Get gene expression from W matrix
pheno.from.W <- function(W, a=0.5, S0=rep(a, nrow(W)), steps=24, measure=4, full=FALSE, loopFUN=internal_loop_R) {
  ans <- loopFUN(W, S0, a, steps, measure)
  if (!full) ans$full <- NULL
  return(ans)
}

#Function to assess after how many steps networks reach "equilibirum"
stablefrom <- function(W, thresh=0.0001, ...) {
  ff <- pheno.from.W(W, ..., full=TRUE)$full
  dif <- rowSums(abs(diff(t(ff))))
  which(dif < thresh)[1]
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
  # ans$eccentricity <- ee$values[2:ncol(M)]/ee$values[1]
  ans$eccentricity <- sqrt(1-(ee$values[2:ncol(M)]/ee$values[1]))
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
  sgn <- if (is.na(target["angle"])) sign(target["cor"]) else sign(target["angle"])
  st <- c(target["size"]/2, target["size"]/2, sgn*target["size"]/10)
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
sim.path <- "../simevolv/bin/Release/Simul_Prog"
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

mutsd.predict <- function(param.multilin, Vm, reps=10, gen=10000, mlim=c(0.001, 0.1), traits=NA, mc.cores=detectCores()-1) {
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

calibrate.multilin <- function(param.multilin, param.wagner, gen=10000, traits=1:2, ...) {
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
oneplot.allellipse <- function(data.dirs, xlim=NULL, ylim=NULL, xlab="Trait A", ylab="Trait B", Sell=TRUE, Mell=TRUE, Gell=FALSE,
                               S.factor=0.00005, M.factor=1,  G.factor=1, xcoord=c(2000, 10000),ycoord=c(-1, 1.6),
                               asp=1, all.gen=FALSE, all.reps=FALSE, legend=TRUE, another_plot=FALSE, ...) {
  
  if (another_plot==TRUE){
  par(fig = c(grconvertX(xcoord, from="user", to="ndc"),
              grconvertY(ycoord, from="user", to="ndc")),
      mar = c(4,6,1,1),
      new = TRUE) 
    }
  
  plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, asp=asp,  ...)

  for (data.dir in data.dirs) {
    
    data.files = list.files(path=data.dir, pattern="simul.*.txt$", full.names=TRUE)
    param.file = list.files(path=data.dir, pattern="*.par",        full.names=TRUE)
    #browser()
    rp <- read.param(param.file)
    if (rp$TYPE_ARCHI %in% c("additive", "multilinear")) {
      mat.col=c(M="darkblue", S="orange", G="darkblue") }
    if (rp$TYPE_ARCHI %in% c("m2")) {
      mat.col=c(M="yellowgreen", S="orange", G="yellowgreen") }
    if (rp$TYPE_ARCHI %in% c("fkl")) {
      mat.col=c(M="maroon2", S="orange", G="maroon2") }
    
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
        G.mat=if(Gell==TRUE) extract.P.matrix(simuls.mean, gen=gen) else NULL,
        G.centre=phen.mean, S.centre=extract.theta(param.file),
        M.factor=M.factor, S.factor=S.factor, G.factor=G.factor,
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
df.fig1 <- function(sims.dir, modulo=pi, all.gen=FALSE){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    print(i)
    files     <- list.files(path = i, full.names=TRUE, pattern = "\\.txt$")
    files <- files[sapply(files, file.size) > 7000]
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
      tt <- fread(ff, data.table=FALSE)
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

df.figsuppG <- function(sims.dir, modulo=pi, all.gen=FALSE){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    print(i)
    files     <- list.files(path = i, full.names=TRUE, pattern = "\\.txt$")
    files <- files[sapply(files, file.size) > 7000]
    param.file.all = list.files(path=i, pattern="\\.par$", full.names=TRUE)
    param.file <- as.character(param.file.all[1])
    stopifnot(length(files) >0)
    
    #read param.par to associate each simu with S properties
    rp <- read.param(param.file)
    mut.rate <- 2*rp$GENET_MUTRATES
    if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    S.mat <- extract.S.matrix(param.file)
    S.ref2 <- matrix.features(S.mat, n.genes=2)[["angle"]][1]
    corrS <- S.mat[1,2]/(sqrt( S.mat[1,1]* S.mat[2,2]))
    mytopos <- lapply(files, function(ff) {
      #print(ff)
      tt <- fread(ff, data.table=FALSE)
      mygens <-rev(if (all.gen==TRUE) tt[,"Gen"] else tt[nrow(tt),"Gen"])
      for (gen in mygens) {
        phen.mean <- extract.P.mean(tt, gen=gen)
        G.mat <- extract.P.matrix(tt, gen=gen)
        M.feat2 <- matrix.features(G.mat,n.genes=2)[["angle"]][1] #M angle
        mean.ma.df <- M.feat2
        M.feat1 <- matrix.features(G.mat,n.genes=2)[["eccentricity"]][1] #M matrix
        corrG <- G.mat[1,2]/(sqrt( G.mat[1,1]* G.mat[2,2]))
        data.gen <- c(ff, gen, S.ref2, M.feat2, M.feat1, corrG, corrS, mean.ma.df )
        filedata <- rbind(filedata, data.gen)
        #create a list of data
      }
      return(as.list(filedata))
    })
    newt <- as.data.frame(rbindlist(mytopos, use.names=FALSE))
    newt[,8] <- mean.angle.pi(as.numeric(newt[,8]))
    setnames(newt, 1:8, c("data.dir", "Gen","ang_S","ang_G", "ecc_G", "corrG","corrS","mean_ang_g"))
    # newt <- ldply(mytopos) # create a data.frame for mytopos
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  simul.df[,2:8] <- lapply( simul.df[,2:8], as.numeric)
  sq_dist <- 1 - (((modulo.all(simul.df$ang_G-simul.df$ang_S, modulo = modulo))^2)/((pi^2)/12))
  simul.df <- as.data.frame(cbind(simul.df, sq_dist)) #df of 8 columns
  return(simul.df)
}


df.topo.wide.m<- function(sims.dir, w_of_6=FALSE, network=FALSE, file_size=100000, all.gen=FALSE){
  simul.df <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    print(i)
    files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
    files <- files[sapply(files, file.size) > file_size]
    param.file = list.files(path=i, pattern="param.par", full.names=TRUE)
    #read param.par to associate each simu with the S "what"
    rp <- read.param(param.file)
    mut.rate <- 2*rp$GENET_MUTRATES
    if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    S.mat <- extract.S.matrix(param.file)
    S.ref2 <- matrix.features(S.mat, n.genes=2)[["angle"]][1] #S matrix
    corrS <- S.mat[1,2]/(sqrt( S.mat[1,1]* S.mat[2,2]))
    mytopos <- lapply(files, function(ff) {
      print(ff)
      tt <- fread(ff, data.table=FALSE) #read.table(ff, header=TRUE)
      mygens <-rev(if (all.gen==TRUE) tt[,"Gen"] else tt[nrow(tt),"Gen"])
      df <- data.frame()
      for (genid in mygens) {
        #M data
        gen <- genid
        M.mat <- extract.M.matrix(tt, gen=genid)
        M.ang_ppi <- matrix.features(M.mat,n.genes=2)[["angle"]][1]+pi #M matrix
        M.ang_mpi <- matrix.features(M.mat,n.genes=2)[["angle"]][1]-pi #M matrix
        M.feat1 <- matrix.features(M.mat,n.genes=2)[["eccentricity"]][1] #M matrix
        M.ang <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
        M.feat3 <- matrix.features(M.mat,n.genes=2)[["size"]][1] #M matrix
        fitness <- extract.fitness(tt, "MFit", gen=genid)
        # opt1 <- extract.fitness(tt, "FitOpt1", gen=gen)[1]
        corrM <- extract.correlation(tt, gen=genid)
        opt2 <- extract.fitness(tt, "FitOpt2", gen=genid)[1]
        if (network==TRUE){
          Wmat <- t(extract.W.matrix(tt)) #W data
          new <- c(ff, gen, M.ang_ppi, S.ref2, M.ang_mpi, M.feat1, M.ang, M.feat3, fitness, corrM, corrS, Wmat )
        }
        if (network!=TRUE){
          new <- c(ff, gen, M.ang_ppi, S.ref2, M.ang_mpi, M.feat1, M.ang, M.feat3, fitness, corrM, corrS )
        }
        df <- rbind(df, new)
      }
      if (network!=TRUE){
        df <- setnames(df[,1:11], c("data.dir","Gen","ang_M_ppi","ang_S","ang_M_mpi","ecc_M","ang_M","siz_M","fitness","corrM","corrS"))
        df[,2:11] <- lapply( df[,2:11], as.numeric)
      }
      if (w_of_6==TRUE & network==TRUE) {
        df <- setnames(df[,1:47], c("data.dir","Gen","ang_M_ppi","ang_S","ang_M_mpi","ecc_M","ang_M","siz_M","fitness","corrM","corrS",
                                                "a_a","b_a","c_a","d_a","e_a","f_a","a_b","b_b","c_b","d_b","e_b","f_b","a_c","b_c","c_c","d_c","e_c","f_c",
                                                "a_d","b_d","c_d","d_d","e_d","f_d","a_e","b_e","c_e","d_e","e_e","f_e","a_f","b_f","c_f","d_f","e_f","f_f"))
        df[,2:47] <- lapply( df[,2:47], as.numeric)
      }
       return(df)
      #create a list of data
    })
    newt <- ldply(mytopos) # create a data.frame for i
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  

  return(simul.df)
}

df.data <- function(sims.dirs, modulo=pi, pattern="simul/", variable="population", network=FALSE, w_of_6=FALSE, file_size=100000, all.gen=FALSE){
  df.datas <- data.frame(NULL)
  for (i in sims.dirs) {
    sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
    df <- df.topo.wide.m(sims.dir, network=network, w_of_6=w_of_6, file_size=file_size, all.gen=all.gen)
    pop <- str_split(i, pattern, n=2, simplify = TRUE)
    df$pop <- sprintf("%s", pop[,2])
    df.datas <- rbind(df.datas, df)
  }
  return(df.datas)
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

