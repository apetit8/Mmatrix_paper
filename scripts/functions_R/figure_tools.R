source("../functions_R/tools.R")

makeTransparent<-function(someColor, alpha=30)
{
        # From https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
        # Author: Nick Sabbe
        # Licence : CC-attribution-SA from the conditions of the website
        newColor<-col2rgb(someColor)
        apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


################################################################################
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

plot.ellipse.dir <- function(data.dir, xlim=NULL, ylim=NULL, xlab="Trait 1", ylab="Trait2", 
                         S.factor=0.00005, M.factor=1, G.factor=1,
                         mat.col=c(G="dodgerblue4", M="darkolivegreen4", S="orange"), 
                         all.gen=FALSE, all.reps=FALSE,  Gell=TRUE, pattern_simu="\\.txt$", ...) {
  data.files = list.files(path=data.dir, pattern=pattern_simu, full.names=TRUE)
  param.file = list.files(path=data.dir, pattern="*.par",        full.names=TRUE)
# stopifnot(length(data.files) >0, length(param.file) == 1)
  stopifnot(length(data.files) >0)
  
  simuls = lapply(data.files, function(x) read.table(x, header=TRUE)) #import list
  simuls.mean = replicate.mean(simuls)
  
  mygens <-rev(if (all.gen) simuls.mean[,"Gen"] else simuls.mean[nrow(simuls.mean),"Gen"])

  for (gen in mygens) {
    if (all.reps) {
      first <- TRUE
      for (simul in simuls) {
        phen.mean <- extract.P.mean(simul, gen=gen)
        draw.ellipses(
          G.mat=extract.P.matrix(simul, gen=gen),
          M.mat=extract.M.matrix(simul, gen=gen),
          S.mat=NULL, G.centre=phen.mean, M.factor=M.factor, S.factor=S.factor, G.factor=G.factor,
          mat.col=makeTransparent(mat.col), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
          add=!(first && gen==mygens[1]), ...)
          first <- FALSE
      }
    }
    # Mean over replicates
    phen.mean <- extract.P.mean(simuls.mean, gen=gen)
    draw.ellipses(
      if (Gell) {
          G.mat=extract.P.matrix(simuls.mean, gen=gen)},
          M.mat=extract.M.matrix(simuls.mean, gen=gen),
          S.mat=if(gen==mygens[1]) extract.S.matrix(param.file[1]) else NULL, 
          G.centre=phen.mean, S.centre=extract.theta(param.file[1]), 
          M.factor=M.factor, S.factor=S.factor, G.factor=G.factor,
          mat.col=mat.col, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
          add=!(!all.reps && gen==mygens[1]), ...)
  }
  if (Gell) {
  legend("topleft", lty=1, col=mat.col, legend=c(paste0("G", if (G.factor != 1) paste0(" (x ", G.factor, ")")), paste0("M", if (M.factor != 1) paste0(" (x ", M.factor, ")")), paste0("S", if (S.factor != 1) paste0(" (x ", S.factor, ")"))))
  }
  else 
    legend("topleft", lty=1, col=mat.col[2:3], legend=c(paste0("M", if (M.factor != 1) paste0(" (x ", M.factor, ")")), paste0("S", if (S.factor != 1) paste0(" (x ", S.factor, ")"))))
}



write.matrix <- function(data.dir) {
  
  data.files = list.files(path=data.dir, pattern="simul.*.txt$", full.names=TRUE)
  param.file = list.files(path=data.dir, pattern="*.par",        full.names=TRUE)
  
  stopifnot(length(data.files) >0, length(param.file) == 1)
  
  simuls = lapply(data.files, function(x) read.table(x, header=TRUE)) #import list
  simuls.mean = replicate.mean(simuls)

  mygens <-simuls.mean[nrow(simuls.mean),"Gen"]
  
  for (gen in mygens) {
    # Phen mean
    phen.mean <- extract.P.mean(simuls.mean, gen=gen)
    # P matrix
    P.mat <- extract.P.matrix(simuls.mean, gen=gen)
    # G matrix
    G.mat <- P.mat # <- diag(nrow(P.mat)) # Only for env variances of 1!!!!
    # M matrix
    M.mat <- extract.M.matrix(simuls.mean, gen=gen)
    #S matrix
    S.mat <- extract.S.matrix(param.file)
    theta <- extract.theta(param.file)
    
    print("M matrix")
    print(M.mat)#M matrix
    print("G matrix")
    print(G.mat)#G matrix
    print("Omega matrix")
    print(S.mat)#S matrix
  
  
}
}

features.extract.matrix <- function(data.dir) {
  
  data.files = list.files(path=data.dir, pattern="simul.*.txt$", full.names=TRUE)
  param.file = list.files(path=data.dir, pattern="*.par",        full.names=TRUE)
  
  stopifnot(length(data.files) >0, length(param.file) == 1)
  
  simuls = lapply(data.files, function(x) read.table(x, header=TRUE)) #import list
  simuls.mean = replicate.mean(simuls)
  
  mygens <-simuls.mean[nrow(simuls.mean),"Gen"]
  
  for (gen in mygens) {
    # Phen mean
    phen.mean <- extract.P.mean(simuls.mean, gen=gen)
    # P matrix
    P.mat <- extract.P.matrix(simuls.mean, gen=gen)
    # G matrix
    G.mat <- P.mat # <- diag(nrow(P.mat)) # Only for env variances of 1!!!!
    # M matrix
    M.mat <- extract.M.matrix(simuls.mean, gen=gen)
    #S matrix
    S.mat <- extract.S.matrix(param.file)
    theta <- extract.theta(param.file)
    
    print("M matrix")
    print(matrix.features(M.mat))#M matrix
    print("G matrix")
    print(matrix.features(G.mat))#G matrix
    print("Omega matrix")
    print(matrix.features(S.mat))#S matrix
    
    plot
  }
}

format.pi.frac <- function(num, den) {
    if (num == 0)
        return(bquote(0))
    if (num == 1)
        if (den == 1)
            return(bquote(pi))
        else
            return(bquote(pi/.(den)))
    else if (num == -1)
        if (den == 1)
            return(bquote(-pi))
        else
            return(bquote(-pi/.(den)))
    else
        if (den == 1)
            return(bquote(.(num)*pi))
        else 
            return(bquote(.(num)*pi/.(den)))
}

plot.features.time <- function(data.dir, what="size", xlim=NULL, ylim=NULL, xlab="Generations", ylab=what, 
                         mat.col=c(G="dodgerblue4", M="darkolivegreen4", S="orange"), 
                         mat.pch=c(G=1, M=1),
                         all.reps=TRUE, ...) {
  stopifnot(what %in% c("size", "eccentricity", "angle", "cor"))
  
  data.files = list.files(path=data.dir, pattern="simul.*.txt$", full.names=TRUE)
  param.file = list.files(path=data.dir, pattern="*.par",        full.names=TRUE)

  stopifnot(length(data.files) >0, length(param.file) == 1)
  
  simuls = lapply(data.files, function(x) read.table(x, header=TRUE)) #import list
  simuls.mean = replicate.mean(simuls)

  mygens <-simuls.mean[,"Gen"]
  rp <- read.param(param.file)
  mut.rate <- 2*rp$GENET_MUTRATES
  if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
  
  G.mean <- sapply(mygens, function(gen) matrix.features(extract.P.matrix(simuls.mean, gen=gen), n.genes=2)[[what]][1])
  M.mean <- sapply(mygens, function(gen) matrix.features(mut.rate*extract.M.matrix(simuls.mean, gen=gen), n.genes=2)[[what]][1])
  S.ref <- matrix.features(extract.S.matrix(param.file), n.genes=2)[[what]][1]

  if (is.null(xlim)) xlim <- range(mygens)
  if (is.null(ylim)) ylim <- if (what=="size") c(S.ref/100000, S.ref) else if (what=="eccentricity") c(0,1) else if (what=="angle") c(S.ref-1.05*pi/2, S.ref+1.05*pi/2) else c(-1,1)

  if (what == "angle") {
    G.mean <- ifelse(G.mean < ylim[1], G.mean + pi, G.mean)
    G.mean <- ifelse(G.mean > ylim[2], G.mean - pi, G.mean)
    M.mean <- ifelse(M.mean < ylim[1], M.mean + pi, M.mean)
    M.mean <- ifelse(M.mean > ylim[2], M.mean - pi, M.mean)
    if (length(S.ref) != 1) browser()
    if(S.ref < ylim[1]) S.ref <- S.ref + pi
    if(S.ref > ylim[2]) S.ref <- S.ref - pi
  }
  
  plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, yaxt=if(what == "angle") "n" else "s", log=if(what=="size") "y" else "", ...)
  if (what=="angle") {
    sq <- seq(-2*pi, 2*pi, by=pi/4)
    sq <- sq[sq > ylim[1] & sq < ylim[2]]
    axis(2, at=sq, label=sapply(sq/pi, function(pif) { rr <- fraction(pif); as.expression(format.pi.frac(rr[1], rr[2])) }))
  }
  
  if (all.reps) {
    for (simul in simuls) {
      # Take the first element of the feature, in case matrix.features returns a vector
      G.points <- sapply(mygens, function(gen) matrix.features(extract.P.matrix(simul, gen=gen), n.genes=2)[[what]][1])
      M.points <- sapply(mygens, function(gen) matrix.features(mut.rate*extract.M.matrix(simul, gen=gen), n.genes=2)[[what]][1])
      if (what == "angle") {
          G.points <- ifelse(G.points < ylim[1], G.points + pi, G.points)
          G.points <- ifelse(G.points > ylim[2], G.points - pi, G.points)
          M.points <- ifelse(M.points < ylim[1], M.points + pi, M.points)
          M.points <- ifelse(M.points > ylim[2], M.points - pi, M.points)
      }
      points(mygens, G.points, pch=mat.pch["G"], col=makeTransparent(mat.col["G"]))
      points(mygens, M.points, pch=mat.pch["M"], col=makeTransparent(mat.col["M"]))
    }
  }

  points(mygens, G.mean, pch=mat.pch["G"], col=mat.col["G"], cex=1.2)
  points(mygens, M.mean, pch=mat.pch["M"], col=mat.col["M"], cex=1.2)
  lines(mygens, rep(S.ref, length(mygens)), col=mat.col["S"])
}


plot.features.onS <- function(data.dirs, what="angle", xlim=NULL, ylim=NULL, xlab="S", ylab="M", 
                              mat.col=c(G="royalblue3", M="darkolivegreen4", S="orange"), 
                              mat.pch=c(G=1, M=16, Gm=17, Mm=17), generation=10000,
                              all.reps=FALSE, Gell=FALSE, ...) {
  stopifnot(what %in% c("size", "eccentricity", "angle", "cor"))
  
  if (is.null(xlim)) xlim <- if (what=="size") c(0,1) else if (what=="eccentricity") c(0,1) else if (what=="angle") c(-pi/2, pi/2) else c(-1,1)
  if (is.null(ylim)) ylim <- if (what=="size") c(0,1) else if (what=="eccentricity") c(0,1) else if (what=="angle") c(-pi/2, pi/2) else c(-1,1)
  plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  
  for (data.dir in data.dirs){
    data.files = list.files(path=data.dir, pattern="simul.*.txt$", full.names=TRUE)
    param.file = list.files(path=data.dir, pattern="*.par",        full.names=TRUE)
    
    stopifnot(length(data.files) >0, length(param.file) == 1)
    
    simuls = lapply(data.files, function(x) read.table(x, header=TRUE)) #import list
    simuls.mean = replicate.mean(simuls)
    
    mygens <-simuls.mean[,"Gen"]
    S.ref <- matrix.features(extract.S.matrix(param.file), n.genes=2)[[what]][1]  
    G.mean <- sapply(S.ref, function(S.ref) matrix.features(extract.P.matrix(simuls.mean, gen=generation), n.genes=2)[[what]][1])
    M.mean <- sapply(S.ref, function(S.ref) matrix.features(extract.M.matrix(simuls.mean, gen=generation), n.genes=2)[[what]][1])
    
    if (all.reps) {
      for (simul in simuls) {
        # Take the first element of the feature, in case matrix.features returns a vector
        if (Gell) {
        G.points <- sapply(S.ref, function(S.ref) matrix.features(extract.P.matrix(simul, gen=generation), n.genes=2)[[what]][1])
        }
        M.points <- sapply(S.ref, function(S.ref) matrix.features(extract.M.matrix(simul, gen=generation), n.genes=2)[[what]][1])
        if (Gell) {
        points(S.ref, G.points, pch=mat.pch["G"], col=makeTransparent(mat.col["G"]))
        }
        points(S.ref, M.points, pch=mat.pch["M"], col=makeTransparent(mat.col["M"]))
      }
    }
    if (Gell) {
    points(S.ref, G.mean, pch=mat.pch["G"], col=mat.col["G"], cex=1.2)
    }
    points(S.ref, M.mean, pch=mat.pch["M"], col=mat.col["M"], cex=1.2)
    #legend("topleft", pch=c(1,16,1), col=mat.col, legend=c("G ", "M ", "S"))
  }
}


plot.features.onS.MvsW <- function(data.dirs1, data.dirs2, what="angle", xlim=NULL, ylim=NULL, xlab="S", ylab=what, 
                              mat.col=c(G="royalblue3", M="darkolivegreen4", S="orange"), 
                              mat.pch=c(G=1, M=1, Gm=6, Mm=6), generation=10000,
                              ...) {
  stopifnot(what %in% c("size", "eccentricity", "angle", "cor"))
  
  if (is.null(xlim)) xlim <- if (what=="size") c(0,1) else if (what=="eccentricity") c(0,1) else if (what=="angle") c(-pi/2, pi/2) else c(-1,1)
  if (is.null(ylim)) ylim <- if (what=="size") c(0,1) else if (what=="eccentricity") c(0,1) else if (what=="angle") c(-pi/2, pi/2) else c(-1,1)
  plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  
  for (data.dir in data.dirs1){
    data.files = list.files(path=data.dir, pattern="simul.*.txt$", full.names=TRUE)
    param.file = list.files(path=data.dir, pattern="*.par",        full.names=TRUE)
    
    stopifnot(length(data.files) >0, length(param.file) == 1)
    
    simuls = lapply(data.files, function(x) read.table(x, header=TRUE)) #import list
    simuls.mean = replicate.mean(simuls)
    
    mygens <-simuls.mean[,"Gen"]
    S.ref <- matrix.features(extract.S.matrix(param.file), n.genes=2)[[what]][1]  
    G.mean <- sapply(S.ref, function(S.ref) matrix.features(extract.P.matrix(simuls.mean, gen=generation), n.genes=2)[[what]][1])
    M.mean <- sapply(S.ref, function(S.ref) matrix.features(extract.M.matrix(simuls.mean, gen=generation), n.genes=2)[[what]][1])
    
    points(S.ref, G.mean, pch=mat.pch["G"], col=mat.col["G"], cex=1.2)
    points(S.ref, M.mean, pch=mat.pch["M"], col=mat.col["M"], cex=1.2)
    
    for (data.dir in data.dirs2){
      data.files = list.files(path=data.dir, pattern="simul.*.txt$", full.names=TRUE)
      param.file = list.files(path=data.dir, pattern="*.par",        full.names=TRUE)
      
      stopifnot(length(data.files) >0, length(param.file) == 1)
      
      simuls = lapply(data.files, function(x) read.table(x, header=TRUE)) #import list
      simuls.mean = replicate.mean(simuls)
      
      mygens <-simuls.mean[,"Gen"]
      S.ref <- matrix.features(extract.S.matrix(param.file), n.genes=2)[[what]][1]  
      G.mean <- sapply(S.ref, function(S.ref) matrix.features(extract.P.matrix(simuls.mean, gen=generation), n.genes=2)[[what]][1])
      M.mean <- sapply(S.ref, function(S.ref) matrix.features(extract.M.matrix(simuls.mean, gen=generation), n.genes=2)[[what]][1])

    points(S.ref, G.mean, pch=mat.pch["Gm"], col=mat.col["G"], cex=1.2)
    points(S.ref, M.mean, pch=mat.pch["Mm"], col=mat.col["M"], cex=1.2)
    legend("topleft", pch=c(16,16, 6, 6, 1), col=c(col=mat.col["G"], col=mat.col["M"], col=mat.col["G"], col=mat.col["M"], col=mat.col["S"]), legend=c("G wagner", "M wagner", "G multi", "M multi","S"))
    }
  }
}


oneplot.allellipse <- function(data.dirs, xlim=NULL, ylim=NULL, xlab="Trait 1", ylab="Trait 2", 
                               S.factor=0.00005, M.factor=1,  G.factor=1, 
                               asp=1, all.gen=FALSE, all.reps=FALSE, ...) {
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
          M.mat=extract.M.matrix(simuls.mean, gen=gen),
          S.mat=if(gen==mygens[1]) extract.S.matrix(param.file) else NULL,
          G.centre=phen.mean, S.centre=extract.theta(param.file),
          M.factor=M.factor, S.factor=S.factor,
          mat.col=mat.col, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
          add=TRUE, ...)
        legend("topleft", lty=1, box.lty=0, bg="transparent", col=c("yellowgreen", "darkblue", "orange"),
               legend=c(paste0("M wagner", if (M.factor != 1) paste0(" (x ", M.factor, ")")), paste0("M multilinear", if (M.factor != 1) paste0(" (x ", M.factor, ")")),paste0("S", if (S.factor != 1) paste0(" (x ", S.factor, ")"))))
    }
  }
}
