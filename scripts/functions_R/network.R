library(Rcpp)
library(ggplot2)
library(plyr)
library(tidyr)
library(rlist)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(data.table)

global.model <- list(
	connect.threshold = 0.01
)


# redefining the base sample function without the interface bug
my.sample <- function(x, size, replace = F, prob = NULL) {
  if (length(x) == 1) return(x)
  base::sample(x, size = size, replace = replace, prob = prob)
}

# redefining the combinat::permn function, for exactly the same reason
my.permn <- function(x, fun=NULL, ...) {
	if (length(x) == 1 ) return(x)
	combinat::permn(x, fun, ...)
}

cppFunction('
	List internal_loop_cpp(const NumericMatrix &W, const NumericVector &S0, double a, unsigned int steps, unsigned int measure) {
		double lambda = (1-a)/a;
		double mu     = 1/(a*(1-a));
		
		NumericMatrix sto (S0.size(), steps+1);
		NumericVector sumx (S0.size());
		NumericVector sumx2 (S0.size()); 

		for (unsigned int i = 0; i < S0.size(); i++)
			sto(i,0) = S0(i);
		for (unsigned int t = 1; t <= steps; t++) {
			for (unsigned int i = 0; i < S0.size(); i++) {
				double tmp = 0.;
				for (unsigned int j = 0; j < S0.size(); j++) {
					tmp += sto(j,t-1) * W(i,j);
				}
				tmp =  1. / (1. + lambda * exp(-mu*tmp));
				
				sto(i,t) = tmp;
				
				if (t > steps-measure) {
					sumx(i) += tmp;
					sumx2(i) += tmp*tmp;
				}
			}
		}
		for (unsigned int i = 0; i < S0.size(); i++) {
			sumx(i) /= static_cast<double>(measure); // sumx(i) now contains the mean
			sumx2(i) /= static_cast<double>(measure);
			sumx2(i) += -sumx(i)*sumx(i); // sumx2(i) now contains the variance
		}
		return List::create(Named("full")=sto, Named("mean")=sumx, Named("var")=sumx2);
	}')


model.M2 <- function(W, S0=rep(a, nrow(W)), a=0.2, steps=24, measure=4, full=FALSE) {
    ans <- internal_loop_cpp(W, S0, a, steps, measure)
	if (!full) ans$full <- NULL
	return(ans)
}

cleanW <- function(W, epsilon=global.model$connect.threshold, ...) {
	# ... are additional arguments to modelM2
	distMat <- matrix(0, ncol=ncol(W), nrow=nrow(W))	
	ref.expr <- model.M2(W, ...)$mean
	
	for (i in 1:nrow(W))
		for (j in 1:ncol(W)) {
			myW <- W
			myW[i,j] <- 0
			new.expr <- model.M2(myW, ...)$mean
			dd <- sqrt(sum((ref.expr-new.expr)^2))
			distMat[i,j] <- dd
		}
	cW <- matrix(0, ncol=ncol(W), nrow=nrow(W))
	cW[distMat > epsilon] <- W[distMat > epsilon]
	cW
}

unique.topo <- function(topo, groups=as.list(1:ncol(topo))) {
	# topo is a n x n square matrix (forced to be a signed matrix)
	# groups is a list of groups of genes, numbered from 1 to the n
	
	# Just in case the user is sloppy, smart guesses:
	topo <- sign(topo)
	if (!all(1:ncol(topo) %in% unlist(groups)))
		groups <- c(groups, as.list((1:ncol(topo))[! 1:ncol(topo) %in% unlist(groups)]))
	
	stopifnot(
		is.matrix(topo), ncol(topo) > 0, ncol(topo) == nrow(topo),
		all(topo %in% c(-1,0,1)),
		length(unlist(groups)) == ncol(topo), all(unlist(groups) %in% 1:ncol(topo))
	)
	
	equiv.topos(topo, groups)[[1]]
}

equiv.topos <- function(topo, groups=as.list(1:ncol(topo)), sorted=TRUE, unique=TRUE) {
	# topo is a n x n square matrix (forced to be a signed matrix)
	# groups is a list of groups of genes, numbered from 1 to the n
	
	# Just in case the user is sloppy, smart guesses:
	topo <- sign(topo)
	if (!all(1:ncol(topo) %in% unlist(groups)))
		groups <- c(groups, as.list((1:ncol(topo))[! 1:ncol(topo) %in% unlist(groups)]))
	
	stopifnot(
		is.matrix(topo), ncol(topo) > 0, ncol(topo) == nrow(topo),
		all(topo %in% c(-1,0,1)),
		length(unlist(groups)) == ncol(topo), all(unlist(groups) %in% 1:ncol(topo))
	)
	
	group.perms <- lapply(groups, function(gr) my.permn(gr))
	glob.perms   <- do.call(expand.grid, lapply(sapply(group.perms, length), seq_len))
	all.perms  <- apply(glob.perms, 1, function(x) unlist(lapply(seq_along(x), function(i) group.perms[[i]][x[i]])))
	
	topo.list <- lapply(as.data.frame(all.perms), function(p) topo[p,p])
	tokeep <- 1:length(topo.list)
	if (unique || sorted) {
		topo.df <- as.data.frame(do.call(rbind, lapply(topo.list, c)))
		if (sorted)
			tokeep <- do.call(order, c(as.list(topo.df[tokeep,]), list(decreasing=TRUE)))
		if (unique)
			tokeep <- tokeep[!duplicated(topo.df[tokeep,])]			
	}
	topo.list[tokeep]
}

plot.network <- function(W, vW=NULL, max=1, col.scale=NULL,
								lwd.arr=2, pos.shift.plus=0.80, pos.shift.minus=0.60, ...) {
	
	circ.arc    <- function(theta1=0, theta2=2*pi, n=100) { tt <- seq(theta1, theta2, length.out=n); cbind(cos(tt), sin(tt)) }
	posit.angle <- function(angle) { angle <- angle %% (2*pi); if (angle > pi/2 && angle <= 3*pi/2) angle <- angle + pi; angle %% (2*pi)}
	mycol       <- function(value) { col.scale[1 + round((value+max)/2/max*(length(col.scale)-1)) ] }
	mylty       <- function(value) { c(0,3,4,2,6,5,1)[1+round(6*value)] }
	
	if (is.null(col.scale))
		col.scale <- colorRampPalette(c("red", "white","black"))(101)
	if (is.null(vW))
		vW <- matrix(1, ncol=ncol(W), nrow=nrow(W))
	
	delta.angle <- 0.4 # angle between two arrows
	arr.dist  <- 0.15 # distance between the group name and the arrows
	self.angle <- 1.4*pi
	ann.text.options <- list(
		pos.shift.plus=pos.shift.plus, 
		pos.shift.minus=pos.shift.minus, 
		text.cex=0.7, 
		col.plus=rev(col.scale)[1], 
		col.minus=col.scale[1], 
		thresh=0.05, 
		digits=2)
	
	par(mar=c(0.1,0.1,2,0.1))
	plot(NULL, xlim=c(-1.2,1.2), ylim=c(-1.2,1.2), axes=FALSE, ann=FALSE, asp=1, ...)
	lg <- ncol(W)
	xy.groups <- cbind(cos(2*pi/lg*(0:(lg-1))), sin(2*pi/lg*(0:(lg-1))))
	
	if (is.null(colnames(W))) colnames(W) <- rownames(W) <- LETTERS[1:ncol(W)] 
	W[W > max] <- max
	W[W < -max] <- -max

	# plots the groups
	text(xy.groups[,1], xy.groups[,2], parse(text=colnames(W)))	

	for (i in 1:lg) {
		for (j in 1:lg) {
			if (i != j) {
				# angle between i and j
				alpha <- atan((xy.groups[j,2]-xy.groups[i,2])/(xy.groups[j,1]-xy.groups[i,1]))
				if (xy.groups[i,1] > xy.groups[j,1]) alpha <- alpha - pi
				
				x.i  <- xy.groups[i,1]+arr.dist*cos(alpha-delta.angle/2)
				y.i  <- xy.groups[i,2]+arr.dist*sin(alpha-delta.angle/2)
				
				x.j  <- xy.groups[j,1]+arr.dist*cos(pi+alpha+delta.angle/2)
				y.j  <- xy.groups[j,2]+arr.dist*sin(pi+alpha+delta.angle/2)
								
				col  <- mycol(W[j,i]) 
				lty  <- mylty(vW[j,i])
				arrows(x0=x.i,  x1=x.j,  y0=y.i,  y1=y.j,  length=0.1,  col=col, lwd=lwd.arr, lty=lty, angle=if(W[j,i] >= 0) 45 else 90)
				
			} else { # i == j
				# angle of i around the circle
				alpha <- (i-1)*2*pi/lg
				# the regulation circle is opposite to the center
				tmp.cc <- circ.arc(alpha-self.angle/2, alpha+self.angle/2)
				cc <- t(t(tmp.cc)*arr.dist+(1+0.5*arr.dist)*xy.groups[i,])
				
				col <- mycol(W[i,i])
				lty  <- mylty(vW[i,i])
				lines(cc[1:(nrow(cc)-1), 1], cc[1:(nrow(cc)-1),2], col=col, lty=lty, lwd=lwd.arr)
				arrows(x0=cc[nrow(cc)-1,1], x1=cc[nrow(cc),1], y0=cc[nrow(cc)-1,2], y1=cc[nrow(cc),2], col=col, length=0.1, lwd=lwd.arr, angle=if(W[j,i] >= 0) 45 else 90)
			}
		}
	}
}

draw.topos <- function(listW, max=1, ...) {
	fun.listmat <- function(lst, FUN) {
		if (length(lst)==1)  {lst <- list(lst[[1]], lst[[1]])}
		n <- length(lst)
		rc <- dim(lst[[1]])
		ar1 <- array(unlist(lst), c(rc, n))
		apply(ar1, c(1, 2), FUN)
	}
	if (is.matrix(listW)) {
		listW <- list(listW)}
	if (is.character(listW)){
		listW <- lapply(strsplit(listW, split="/"), function(ss) matrix(as.numeric(ss), ncol=sqrt(length(ss))))}
	
	listW <- lapply(listW, function(W) { W[W < -max] <- -max; W[W > max] <- max; W })
	mm <- fun.listmat(listW, mean)
	vv <- fun.listmat(listW, var)
	mxvar <- max^2
	plot.network(W=mm, vW=1-vv/mxvar, max=max, ...)
}

#create the dataframe used to analyse topo
df.feat.and.topo <- function(sims.dir, of, what, where, threshold, c){
  simul.topo1 <- data.frame()
  
  #collect all the data
  for (i in sims.dir) {
    files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
    param.file = list.files(path=i, pattern="param.par", full.names=TRUE)
    stopifnot(length(files) >0, length(param.file) == 1)
    
    #read param.par to associate each simu with the S "what"
    rp <- read.param(param.file)
    mut.rate <- 2*rp$GENET_MUTRATES
    if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    if (what == "ecc") {
      S.ref <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["eccentricity"]][1] 
    }
    else { 
      S.ref <- matrix.features(extract.S.matrix(param.file), n.genes=2)[[what]][1]
    }

    
    mytopos <- lapply(files, function(ff) {
      #print(ff)
      tt <- read.table(ff, header=TRUE)
      Wmat <- extract.W.matrix(tt)
      cWmat <- cleanW(Wmat, a=c, epsilon=threshold) #use the threshold to group the topology who are the same 
      untopo <-unique.topo(cWmat, groups=list(1:2,3:4)) #group the topology who are exancheable
      
      #"what" of G and M
      gen <-tt[nrow(tt),"Gen"]
      P.mat <- extract.P.matrix(tt, gen=gen)
      M.mat <- extract.M.matrix(tt, gen=gen)
      if (what == "ecc") {
        M.feat <- matrix.features(M.mat,n.genes=2)[["eccentricity"]][1] #M matrix
        G.feat <- matrix.features(P.mat,n.genes=2)[["eccentricity"]][1]#G matrix   
      }
      else { 
        M.feat <- matrix.features(M.mat,n.genes=2)[[what]][1] #M matrix
        G.feat <- matrix.features(P.mat,n.genes=2)[[what]][1]#G matrix
      }
       new <- c(ff, S.ref, M.feat, G.feat, untopo) #create a list of wich topology goes with wich simulation file
    })
    newt <- ldply(mytopos) # create a data.frame for the i "what"
    simul.topo1 <- rbind(simul.topo1, newt) #paste the data of the i "what" in a data.frame with the others
    }
  
  names(simul.topo1)[names(simul.topo1) == "V1"] <- "data.dir"
  names(simul.topo1)[names(simul.topo1) == "V2"] <- "what_S"
  names(simul.topo1)[names(simul.topo1) == "V3"] <- "what_M"
  names(simul.topo1)[names(simul.topo1) == "V4"] <- "what_G"
  
  #concatenate  the topology in 1 column
  simul.topo1$topo = paste(simul.topo1[,5],simul.topo1[,6],simul.topo1[,7],simul.topo1[,8],simul.topo1[,9],simul.topo1[,10],simul.topo1[,11],simul.topo1[,12],simul.topo1[,13],simul.topo1[,14],simul.topo1[,15],simul.topo1[,16],simul.topo1[,17],simul.topo1[,18],simul.topo1[,19],simul.topo1[,20])

  #create a colomn to link angles with topology
  library(data.table) ## v 1.9.6+ 
  setDT(simul.topo1)[, paste0("name", 1:2) := tstrsplit(data.dir, sprintf("/%s_", what))]
  setDT(simul.topo1)[, paste0("name", 1:2) := tstrsplit(tolower(name2), sprintf(".par/simul%s", what))]
  names(simul.topo1)[names(simul.topo1) == "name1"] <- "Value"
  names(simul.topo1)[names(simul.topo1) == "name2"] <- "Simul"
  
  simul.topo1$Value <- as.numeric(as.character(simul.topo1$Value))
  simul.topo1$what_S <- as.numeric(as.character(simul.topo1$what_S))
  simul.topo1$what_G <- as.numeric(as.character(simul.topo1$what_G))
  simul.topo1$what_M <- as.numeric(as.character(simul.topo1$what_M))
  
  # complexity1 - number of connections
  comp <- as.data.frame( lapply((simul.topo1[,5:20]), (as.numeric) ))
  simul.topo1$complexity1 <- apply(comp, 1, function(x) sum(x %in% c(1,-1)))
  # complexity2 - number of connections with non selected traits
  comp <- as.data.frame( lapply((simul.topo1 %>% select(8,9,12:20)), (as.numeric) ))
  simul.topo1$complexity2 <- apply(comp, 1, function(x) sum(x %in% c(1,-1)))
  #  complexity3 - 1 count
  comp <- as.data.frame( lapply((simul.topo1[,5:20]), (as.numeric) ))
  simul.topo1$complexity3 <- apply(comp, 1, function(x) sum(x %in% c(1)))
  # complexity4 - -1 count
  comp <- as.data.frame( lapply((simul.topo1[,5:20]), (as.numeric) ))
  simul.topo1$complexity4 <- apply(comp, 1, function(x) sum(x %in% c(-1)))

  #order topologies
  index   <- simul.topo1$complexity1
  to.sort <- as.factor(simul.topo1$topo)
  order(index)
  simul.topo1$topo <- factor(simul.topo1$topo, levels=unique(to.sort[order(index)]))
  
  #distances
  simul.topo1$dist_GS <- simul.topo1$what_G-simul.topo1$what_S
  simul.topo1$dist_MS <- simul.topo1$what_M-simul.topo1$what_S
  simul.topo1$absdist_GS <- abs(modulopi(simul.topo1$what_S-simul.topo1$what_G))
  simul.topo1$absdist_MS <- abs(modulopi(simul.topo1$what_S-simul.topo1$what_M))
  
  #add frequency column
  freq <- data.frame()
  for (i in simul.topo1$topo) {
    occ <- simul.topo1[which(simul.topo1$topo == i), ]
    freq <- rbind(freq, nrow(occ))
  }
  simul.topo1$Frequency <- freq
  #use this to create new column with the important topologies and "others"
  matter <- data.frame()
  for (i in c(1:nrow(simul.topo1))) {
    if (simul.topo1[i,32] >= nrow(simul.topo1)*0.020) {
      newt <- c(simul.topo1[i,21])
    }
    else {
      newt <- c("others")
    }
    matter <- rbind(unname(matter), unname(newt))
  }
  simul.topo1$matter <- matter

  return(simul.topo1)
  
}

#draw mean topo for angles in direction
draw.topo.angle.from.dir <- function(sims.dir, threshold){
  for (i in sims.dir) {
    files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
    mytopos <- lapply(files, function(ff) {
      #print(ff)
      tt <- read.table(ff, header=TRUE)
      Wmat <- extract.W.matrix(tt)
      cWmat <- cleanW(Wmat, a=c, epsilon=threshold) #use the threshold to group the topology who are the same 
      untopo <-unique.topo(cWmat) #group the topology who are exancheable
    })
    draw.topos(mytopos)
  }
}


df.topo.raw<- function(sims.dir, w_of_4=TRUE, network=TRUE){
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
    S.ref1 <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["eccentricity"]][1]
    S.ref2 <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["angle"]][1]
    S.ref3 <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["size"]][1]
    
    mytopos <- lapply(files, function(ff) {
      #print(ff)
      tt <- read.table(ff, header=TRUE)
      #M data
      gen <-tt[nrow(tt),"Gen"]
      M.mat <- extract.M.matrix(tt, gen=gen)
      M.feat1 <- matrix.features(M.mat,n.genes=2)[["eccentricity"]][1] #M matrix
      M.feat2 <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
      M.feat3 <- matrix.features(M.mat,n.genes=2)[["size"]][1] #M matrix
      mean.ma.df <- M.feat2
      if (network==TRUE){
        Wmat <- extract.W.matrix(tt) #W data
        new <- c(ff, S.ref1, S.ref2, S.ref3, M.feat1, M.feat2, M.feat3, mean.ma.df, Wmat )
      }
      if (network!=TRUE){
        new <- c(ff, S.ref1, S.ref2, S.ref3, M.feat1, M.feat2, M.feat3, mean.ma.df )
      }
      return(new)
      #create a list of data
    })
    newt <- ldply(mytopos) # create a data.frame for i
    newt[,8] <- mean.angle.pi(as.numeric(newt[,8]))
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  simul.df[,2:8] <- lapply( simul.df[,2:8], as.numeric)
  
  names(simul.df)[names(simul.df) == "V1"] <- "data.dir"
  names(simul.df)[names(simul.df) == "V2"] <- "ecc_S"
  names(simul.df)[names(simul.df) == "V3"] <- "ang_S"
  names(simul.df)[names(simul.df) == "V4"] <- "siz_S"
  names(simul.df)[names(simul.df) == "V5"] <- "ecc_M"
  names(simul.df)[names(simul.df) == "V6"] <- "ang_M"
  names(simul.df)[names(simul.df) == "V7"] <- "siz_M"
  names(simul.df)[names(simul.df) == "V8"] <- "mean_ang_M"
  
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

df.raster.map <- function(sims.dir, modulo=pi, all.gen=FALSE){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    files     <- list.files(path = i, full.names=TRUE, pattern = "\\.txt$")
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
  sq_dist <- 1 - ((modulo.all(simul.df$ang_M-simul.df$ang_S, modulo = modulo))^2)/((pi^2)/12)
  simul.df <- as.data.frame(cbind(simul.df, sq_dist))
  return(simul.df)
  
}
