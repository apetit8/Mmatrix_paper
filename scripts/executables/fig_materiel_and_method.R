#Schematic plots to explain the relationship between e(M), alpha(M), and corr(M)

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

precision <- 101

rr <- seq(-0.99, 0.99, length.out=precision)
aa <- seq(-pi/2, pi/2, length.out=precision)
ee <- seq(0, 1, length.out=precision)

r.Fae <- outer(aa, ee, function(aaa, eee) mapply(aaa, eee, FUN=function(a, e) matrix.features(matrix2.from.features(angle=a, eccentricity=e, size=1))$cor))
a.Fer <- outer(ee, rr, function(eee, rrr) mapply(eee, rrr, FUN=function(e, r) matrix.features(matrix2.from.features(eccentricity=e, cor=r, size=1))$angle[1]))
e.Far <- outer(aa, rr, function(aaa, rrr) mapply(aaa, rrr, FUN=function(a, r) matrix.features(matrix2.from.features(angle=a, cor=r, size=1))$eccentricity))

col.r <- hcl.colors(1024, "PRGn", rev = TRUE)
col.e <- hcl.colors(1024, "Heat2", rev=TRUE)

# 
# pdf("ellipsgeom.pdf", width=12, height=4)
# layout(t(1:3))
# 
# image(aa, ee, r.Fae, xlab="Angle", ylab="Eccentricity", main="Correlation", col=col.r)
# contour(aa, ee, r.Fae, add=TRUE)
# 
# image(ee, rr, a.Fer, xlab="Eccentricity", ylab="Correlation", main="Angle", col=col.r)
# contour(ee, rr, a.Fer, add=TRUE, levels=seq(-3*pi/4, 3*pi/4, length=5))
# 
# image(aa, rr, e.Far, xlab="Angle", ylab="Correlation", main="Eccentricity", col=col.e)
# contour(aa, rr, e.Far, add=TRUE)
# dev.off()

plot.new()
pdf("figures/fig_matmet_curvecor.pdf", width=5, height=5)
plot(NULL, xlim=range(aa), ylim=c(-1,1), xlab=expression(alpha(M)), ylab="Correlation r(M)", xaxt="n", cex.axis=1.2, cex.lab=1.2)
e.test <- c(0.2, 0.7, 0.85, 0.95)
for (ie in seq_along(e.test)) {
  lines(aa, sapply(aa, function(a) matrix.features(matrix2.from.features(angle=a, eccentricity=e.test[ie], size=1))$cor), col=ie)
}
axis(side=1, at=c(-pi, -pi/2, -pi/4, 0 , pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi), mgp = c(1.75, 0.75, 0), cex.axis=1.2)
legend("bottomright", lty=1, col=seq_along(e.test), legend=paste("e = ", e.test))
dev.off()







#

##############
ang1 <- -pi/4
ang2 <- pi/4
ang3 <- 0
ang3.1 <- -pi/2
ang3.2 <- pi/2
ang4 <- pi/8
##############

alpha1 = function(x){ang1}
alpha2 = function(x){ang2}
alpha3 = function(x){ang3}
alpha3.1 = function(x){ang3.1}
alpha3.2 = function(x){ang3.2}
alpha4 = function(x){ang4}

g1 <- ggplot(data.frame(x=c(0, 1)), aes(x=x)) + 
  stat_function(fun=alpha1, colour="blue", show.legend = FALSE)+
  stat_function(fun=alpha2, colour="red", show.legend = FALSE)+
  stat_function(fun=alpha3, colour="#000000", show.legend = FALSE)+
  stat_function(fun=alpha3.1, colour="#000000", show.legend = FALSE)+
  stat_function(fun=alpha3.2, colour="#000000", show.legend = FALSE)+
  stat_function(fun=alpha4, colour="darkgreen", show.legend = FALSE)+
  scale_y_continuous(limits = c(-pi/2-0.1,pi/2+0.1), breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  labs(y=expression(paste(alpha, "(M)")), x=expression(paste("e(M)")))+
  theme_bw() + ggtitle("A")
#
#
corr1 = function(x){-x}
corr2 = function(x){x}
corr3 = function(x){0*x}
corr4 = function(x){x*0.5}

g2 <- ggplot(data.frame(x=c(0, 1)), aes(x=x)) + 
  stat_function(fun=corr1, colour="blue", show.legend = FALSE)+
  stat_function(fun=corr2, colour="red", show.legend = FALSE)+
  stat_function(fun=corr3, colour="#000000", show.legend = FALSE)+
  stat_function(fun=corr4, colour="darkgreen", show.legend = FALSE)+
  scale_y_continuous(limits = c(-1,1))+
  labs(y=expression(paste("r(M)")), x=expression(paste("e(M)")))+theme_bw() + ggtitle(" ")

################################################################################


ecc1 = function(x){sqrt(1-0)}
ecc2 = function(x){sqrt(1-0.25)}
ecc3 = function(x){sqrt(1-0.99)}
ecc4 = function(x){sqrt(1-0.5)}

g3 <- ggplot(data.frame(x=c(-pi/2, pi/2)), aes(x=x)) + 
  stat_function(fun=ecc1, colour="blue", show.legend = FALSE)+
  stat_function(fun=ecc2, colour="red", show.legend = FALSE)+
  stat_function(fun=ecc3, colour="#000000", show.legend = FALSE)+
  stat_function(fun=ecc4, colour="darkgreen", show.legend = FALSE)+
  scale_x_continuous(limits = c(-pi/2-0.1,pi/2+0.1), breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  labs(x=expression(paste(alpha, "(M)")), y=expression(paste("e(M)")))+theme_bw() + ggtitle("B")




alpha1 = function(x){1*sin(x = (pi + 2*x))}
alpha2 = function(x){0.5*sin(x = (pi + 2*x))}
alpha3 = function(x){0}
alpha4 = function(x){0.25*sin(x = (pi + 2*x))}

g4 <- ggplot(data.frame(x=c(-pi/2, pi/2)), aes(x=x)) + 
  stat_function(fun=alpha1, colour="blue", show.legend = FALSE)+
  stat_function(fun=alpha2, colour="red", show.legend = FALSE)+
  stat_function(fun=alpha3, colour="#000000", show.legend = FALSE)+
  stat_function(fun=alpha4, colour="darkgreen", show.legend = FALSE)+
  scale_x_continuous(limits = c(-pi/2-0.1,pi/2+0.1), breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  labs(x=expression(paste(alpha, "(M)")), y=expression(paste("r(M)")))+theme_bw() + ggtitle(" ")


cairo_pdf("../../figures/Mmatrix_schema.pdf", width=8, height=8)
grid.arrange(
  g1,g2,g3,g4,
  ncol = 2,
  nrow = 2,
  widths = c(1,1),
  clip = FALSE
)
dev.off()






