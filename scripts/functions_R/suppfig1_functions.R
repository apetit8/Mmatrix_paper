source("../functions_R/All_functions.R")
library(optimx)
library(parallel)
mc.cores <- min(detectCores()-1, 60)

########Tools###################################################################

make.sym.mat <- function(m) {
	# upper tri copied in lower tri
	m[lower.tri(m)] = t(m)[lower.tri(m)]
	m
}

M.matrix <- function(W, fixed=NULL, sigmaM=0.1, ...) {
	
	ij <- seq_along(W)
	Ps <- mclapply(ij[!ij %in% fixed], function(i) {
			myWm <- myWp <- W
			myWm[i] <- W[i] - sigmaM
			myWp[i] <- W[i] + sigmaM
			return(rbind(model.M2(myWm, ...)$mean, model.M2(myWp, ...)$mean))
	}, mc.cores=1)
	var(do.call(rbind, Ps))
}


S.matrix <- function(s, optalpha, eccentricity=0.7) {
	.matrix2.optim <- function(target) {
		stopifnot(all(names(target) %in% c("cor","angle","size","eccentricity")))
		stopifnot(length(target) == 3)
		
		.matrix2.targetdist <- function(pp, target) {
			if (pp["v1"] <= 0 || pp["v2"] <= 0) return(Inf)
			mf <- matrix.features(matrix(c(pp["v1"], pp["c"], pp["c"], pp["v2"]), ncol=2),n.genes=2)
			mf <- sapply(mf, function(x) x[1])
			sum(sapply(names(target), function(nn) mf[[nn]][1]-target[nn])^2)
		}
  
		st <- c(target["size"]/2, target["size"]/2, sign(target["angle"])*target["size"]/10)
		names(st) <- c("v1","v2","c")
		ans <- optimx::optimx(par=st, .matrix2.targetdist, method="Nelder-Mead", target=target, itnmax=2000)
		if (ans$convcode != 0) 
			warning("Convergence failed for ", paste0(names(target), "=", target, collapse="  "))
		matrix(c(ans$v1, ans$c, ans$c, ans$v2), ncol=2)
	}
	.matrix2.optim(c(size=sum(1/s[1:2]), angle=optalpha, eccentricity=eccentricity))
}

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
########Simevolv################################################################

sigma.M2 <- function(x, a) {
  1. / (1. + exp((-x/(a*(a-1)))+log(1/a-1)))
}

sigma.M2p <- function(x, lambda, mu) {
  1. / (1. + lambda * exp(-mu*x))
}

suppressMessages(library(compiler))
sigma.M2c <- cmpfun(sigma.M2p)

#~ library(Rcpp)

#~ cppFunction('
#~ NumericVector sigma_M2p(NumericVector x, double aam1, double l1am1) {
#~     NumericVector ans(x.size());
#~     // double aam1 = a*(1.-a);
#~     // double l1am1 = log(1./a-1.);
#~     for (int i = 0; i < x.size(); i++) {
#~         ans[i] = 1. / (1. + exp((-x[i]/aam1)+l1am1));
#~     }
#~     return ans;
#~ }
#~ ')

#~ cppFunction('
#~ NumericVector sigma_M2(NumericVector x, double a) {
#~     NumericVector ans(x.size());
#~     double aam1 = a*(1.-a);
#~     double l1am1 = log(1./a-1.);
#~     for (int i = 0; i < x.size(); i++) {
#~         ans[i] = 1. / (1. + exp((-x[i]/aam1)+l1am1));
#~     }
#~     return ans;
#~ }
#~ ')

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

library(Rcpp)
library(inline, quietly=TRUE)

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


model.M2 <- function(W, a=0.5, S0=rep(a, nrow(W)), steps=20, measure=4, full=FALSE, loopFUN=internal_loop_cpp) {
  ans <- loopFUN(W, S0, a, steps, measure)
  if (!full) ans$full <- NULL
  return(ans)
}

model.M2 <- cmpfun(model.M2)


########Optimization pipeline###################################################
n.genes <- 6
optim.model <- c("alpha", "fitness")[1]
#~ optim.model <- "fitness"

default.set <- list(
  s = c(10, 10, 10, 10, 0, 0),
  sp = rep(100, n.genes),
  sa = 1, 
  se = 1, 
  ss = 100,
  optP = c(0.5, 0.5, 0.7, 0.7, 0.5, 0.5),
  optsize = 0.026
)

default.par <- list(
  a=0.5, 
  steps=20, 
  measure=4
)

propW <- function(W, ...) {
  phen <- model.M2(W, ...)
  mut  <- M.matrix(W, ...)
  mut.prop <- matrix.features(mut, n.genes=2)
  list(P=phen$mean, varP=phen$var, M=mut, alphaM=mut.prop$angle[1], excent=1-mut.prop$eccentricity[1], size=mut.prop$size)
}

minus.log.fitness.P <- function(P, varP, alphaM, excent, size, s, sp, sa, se, ss, optP, optalpha, optsize) {
  ff <- 
    sum(s*(P-optP)^2) + sum(sp*varP) + sa*(modulopi(alphaM-optalpha)^2) + se*(1-excent)^2 + ss*(size-optsize)^2
}

minus.log.fitness.M <- function(P, varP, optP, s, sp, M, S) {	
  M <- M[1:2,1:2]
  sum(s*(P-optP)^2) + sum(sp*varP) - sqrt(det(S %*% solve(S+M) %*% M)/det(M))
}

minus.log.fitness.W <- function(W, s, sp, sa, se, ss, optP, optalpha, optsize, S=NULL, ...) {
  pW <- propW(W, ...)
  if (optim.model == "alpha") {
    minus.log.fitness.P(P=pW$P, varP=pW$varP, alphaM=pW$alphaM, excent=pW$excent, size=pW$size, s=s, sp=sp, sa=sa, se=se, ss=ss, optP=optP, optalpha=optalpha, optsize=optsize)
  } else {
    if (is.null(S)) 
      S <- S.matrix(s, optalpha)
    minus.log.fitness.M(P=pW$P, varP=pW$varP, optP=optP, s=s, sp=sp, M=pW$M, S=S)
  }
}

optimW <- function(optalpha, templateW, sigP0=0) {
  ng <- ncol(templateW)
  pp0 <- setNames(rnorm(ng*ng, mean=0, sd=sigP0), nm=outer(seq_len(ng), seq_len(ng), function(w1, w2) paste0("W", w1, w2)))
  pp0 <- pp0[is.na(templateW)]
  
  myS <- S.matrix(default.set$s, optalpha)
  
  fn <- function(par) {
    myW <- templateW
    myW[is.na(templateW)] <- par
    ans <- 
      minus.log.fitness.W(myW, 
                          s=default.set$s, sp=default.set$sp, sa=default.set$sa, se=default.set$se, ss=default.set$ss, 
                          optP=default.set$optP, optalpha=optalpha, optsize=default.set$optsize,
                          a=default.par$a, steps=default.par$steps, measure=default.par$measure, S=myS)
    #~ 		mlw.hist <<- c(mlw.hist, propW(myW)$alphaM)
    ans
  }
  oo <- optimx(pp0, fn, method="L-BFGS-B", lower=rep(-1, length(pp0)), upper=rep(1, length(pp0)), itnmax=500)
  ansW <- templateW
  ansW[is.na(templateW)] <- unlist(oo[seq_len(length(pp0))])
  attr(ansW, "convcode") <- oo$convcode
  ansW
}

optimW.alpharange <- function(templateW, sigP0=0, n.points=101) {
  alphas <- seq(-pi/2, pi/2, length.out=n.points)
  ll <- mclapply(alphas, optimW, templateW=templateW, sigP0=sigP0, mc.cores=mc.cores)
  pp <- lapply(ll, function(W) propW(W, a=default.par$a, steps=default.par$steps, measure=default.par$measure))
  list(
    alphaS = alphas,
    alphaM = sapply(pp, "[[", "alphaM"),
    excent = sapply(pp, "[[", "excent"),
    size   = sapply(pp, "[[", "size"),
    W      = ll
  )
}

########Figs####################################################################

.plot_oa <- function(oa) {
  layout(t(1:3))
  plot(oa$alphaS, oa$alphaM, xlab="alpha S", ylab="alpha M")
  plot(oa$alphaS, oa$excent, xlab="alpha S", ylab="eccentricity", ylim=c(0,1))
  plot(oa$alphaS, oa$size, xlab="alpha S", ylab="size", ylim=c(0, max(oa$size)))
}

.make_df <- function(oa, label) {
  ans <- data.frame(series=label, alphaS=oa$alphaS, alphaM=oa$alphaM, eccentricity=oa$excent, size=oa$size)
  wws <- do.call(rbind, lapply(oa$W, function(w) c(t(w))))
  colnames(wws) <- paste0("W", c(t(outer(seq_len(ncol(oa$W[[1]])), seq_len(nrow(oa$W[[1]])), paste0))))
  cbind(ans, wws)
}