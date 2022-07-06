source("../functions_R/All_functions.R")
#########################################
sims.dirs <- list.dirs("../../simul/fig_1abc", recursive = FALSE)
angle = c(0.785, 0, -0.393)
of        <- "fig1"
modulo <- pi
#####################

M.factor <- 1500
G.factor <- 0.25
S.factor <- 1

#Df with generations
df.m.s.gen <- data.frame()
for (i in angle) {
  sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle%s.parangle", i), full.names=TRUE)
  df <- df.fig1(sims.dir, all.gen=TRUE)
  pop <- str_split(df$data.dir, "../../simul/fig_1abc/", n=2, simplify = TRUE)
  df[,9] <- sprintf("%s", pop[,2])
  model <- str_split(pop[,2], "/simuangle", n=2, simplify = TRUE)
  df[,10] <- sprintf("%s", model[,1])
  df.m.s.gen <- rbind(df.m.s.gen, df) 
}
df.m <- subset(df.m.s.gen, V10 == "m") #Df of Multilinear simulations
df.w <- subset(df.m.s.gen, V10 == "w") #Df of Wagner simulations

colors <- c("darkblue", "yellowgreen")


# png(file="../../figures/fig1_part1.png", width=400, height=400)
cairo_pdf("../../figures/fig1_a.pdf", width=6, height=6)
  sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle%s.parangle", 0.785), full.names=TRUE)
  mar=c(0,0,0,0)
  dfang1 <- subset(df.m.s.gen, round(ang_S, 1) == round(0.785, 1))
  dfm <- subset(df.m, round(ang_S, 1) == round(0.785, 1))
  dfw <- subset(df.w, round(ang_S, 1) == round(0.785, 1))
  plot(dfang1$Gen, dfang1$ang_M, ylim =c(-pi/2,pi/2), xlim = c(min(dfang1$Gen), max(dfang1$Gen)),
       main="", yaxt="n", ylab = expression(paste("Mutational effects direction (",alpha, "M)")), xlab = "Generation", col=alpha(colors[factor(dfang1$V10)],0.2),mar=c(1, 1, 0, 1), mgp = c(1.75, 0.75, 0))
  # for (i in unique(as.factor(dfang1$data.dir))) {
  #   lines(  subset(dfang1, data.dir==i)$Gen,  subset(dfang1, data.dir==i)$ang_M, col=alpha(colors[factor(subset(dfang1, data.dir==i)$V10)],0.2)  )
  # }
  df500 <- subset(dfang1, Gen >= 500)
  bymodel <- by(df500$ang_M, list(df500$Gen, df500$V10), FUN=mean.angle.pi)
  lines(as.numeric(rownames(bymodel)), bymodel[,"m"], col="darkblue", lwd = 3)
  lines(as.numeric(rownames(bymodel)), bymodel[,"w"], col="yellowgreen", lwd = 3)
  lines(dfang1$Gen, dfang1$ang_S, col="orange", type = "l", lty=3, lwd = 2) #S orientation
  axis(side=2, at=c(-pi, -pi/2, -pi/4, 0 , pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi), mgp = c(1.75, 0.75, 0))
  legend("bottomleft", lty=1, box.lty=0,  bg="transparent", col=c("yellowgreen", "darkblue", "orange"),
         legend=c(paste0("M GRN"), paste0("M multilinear"),paste0("S")))
  
  
  oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, another_plot=TRUE,
                     xlim=c(-7,7), ylim=c(-7,7), all.reps=FALSE, xlab="", ylab="",xcoord=c(1000, 10000),ycoord=c(-2.5, 0.5),
                     yaxt = "n", xaxt = "n", legend=legend, asp=1, axes = FALSE,mgp = c(0, 0, 0))
dev.off()
########

# png(file="../../figures/fig1_part2.png", width=400, height=400) 
cairo_pdf("../../figures/fig1_b.pdf", width=6, height=6)
  sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle%s.parangle", 0), full.names=TRUE)
  dfang1 <- subset(df.m.s.gen, round(ang_S, 1) == round(0, 1))
  dfm <- subset(df.m, round(ang_S, 1) == round(0, 1))
  dfw <- subset(df.w, round(ang_S, 1) == round(0, 1))
  plot(dfang1$Gen, dfang1$ang_M, ylim =c(-pi/2,pi/2), xlim = c(min(dfang1$Gen), max(dfang1$Gen)),
       main="", yaxt="n", ylab = expression(paste(alpha, "M")), xlab = "Generation", col=alpha(colors[factor(dfang1$V10)],0.2), mgp = c(1.75, 0.75, 0))
  df500 <- subset(dfang1, Gen >= 500)
  bymodel <- by(df500$ang_M, list(df500$Gen, df500$V10), FUN=mean.angle.pi)
  lines(as.numeric(rownames(bymodel)), bymodel[,"m"], col="darkblue", lwd = 3)
  lines(as.numeric(rownames(bymodel)), bymodel[,"w"], col="yellowgreen", lwd = 3)
  lines(dfang1$Gen, dfang1$ang_S, col="orange", type = "l", lty=3, lwd = 2) #S orientation
  axis(side=2, at=c(-pi, -pi/2, -pi/4, 0 , pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi), mgp = c(1.75, 0.75, 0))

  oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, another_plot=TRUE,
                     xlim=c(-7,7), ylim=c(-7,7), all.reps=FALSE, xlab="", ylab="",xcoord=c(0, 10000),ycoord=c(-2.4, 0.3),
                     yaxt = "n", xaxt = "n", legend=legend, asp=1, axes = FALSE,mgp = c(0, 0, 0))
dev.off()
#######

# png(file="../../figures/fig1_part3.png", width=400, height=400)
cairo_pdf("../../figures/fig1_c.pdf", width=6, height=6)
  sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle%s.parangle", -0.393), full.names=TRUE)
  dfang1 <- subset(df.m.s.gen, round(ang_S, 1) == round(-0.393, 1))
  dfm <- subset(df.m, round(ang_S, 1) == round(-0.393, 1))
  dfw <- subset(df.w, round(ang_S, 1) == round(-0.393, 1))
  plot(dfang1$Gen, dfang1$ang_M, ylim =c(-pi/2,pi/2), xlim = c(min(dfang1$Gen), max(dfang1$Gen)),
       main=NULL, yaxt="n", ylab = expression(paste(alpha, "M")), xlab = "Generation", col=alpha(colors[factor(dfang1$V10)],0.2), mgp = c(1.75, 0.75, 0))
  df500 <- subset(dfang1, Gen >= 500)
  bymodel <- by(df500$ang_M, list(df500$Gen, df500$V10), FUN=mean.angle.pi)
  lines(as.numeric(rownames(bymodel)), bymodel[,"m"], col="darkblue", lwd = 3)
  lines(as.numeric(rownames(bymodel)), bymodel[,"w"], col="yellowgreen", lwd = 3)
  lines(dfang1$Gen, dfang1$ang_S, col="orange", type = "l", lty=3, lwd = 2) #S orientation
  axis(side=2, at=c(-pi, -pi/2, -pi/4, 0 , pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi), mgp = c(1.75, 0.75, 0))

  oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, another_plot=TRUE,
                     xlim=c(-7,7), ylim=c(-7,7), all.reps=FALSE, xlab="", ylab="", xcoord=c(0, 10000),ycoord=c(-1, 2),
                     yaxt = "n", xaxt = "n", legend=legend, asp=1, axes = FALSE,mgp = c(0, 0, 0))

dev.off()
# dev.off()
# print("Ellipses done !")


# 
# cairo_pdf("../../figures/fig1.pdf", width=6.5, height=8)
# laymat <- matrix(c(1, 3, 5, 2, 4, 6), nrow = 3, ncol = 2)
# layout(mat = laymat) # Widths of the two columns
# legend <- TRUE
# for (i in angle) {
#   sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle%s.parangle", i), full.names=TRUE)
#   par(mar = c(3.2, 3 , 0.2, 1))
#   
#   oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, 
#                      # main= paste0("\u03BE", "\u03B1", sprintf(" : m = %s, w = %s", mse.m, mse.w)),
#                      xlim=c(-7,7), ylim=c(-7,7), all.reps=FALSE,
#                      yaxt = "n", xaxt = "n", mgp = c(1, 1, 0), legend=legend, asp=1)
#   dfang1 <- subset(df.m.s.gen, round(ang_S, 1) == round(i, 1))
#   dfm <- subset(df.m, round(ang_S, 1) == round(i, 1))
#   dfw <- subset(df.w, round(ang_S, 1) == round(i, 1))
#   plot(dfang1$Gen, dfang1$ang_M, ylim =c(-pi/2,pi/2), xlim = c(min(dfang1$Gen), max(dfang1$Gen)),
#        main="", yaxt="n", ylab = "M direction", xlab = "Generations", col=alpha(colors[factor(dfang1$V10)],0.2), mgp = c(1.75, 0.75, 0))
#   # lines(lowess(dfm$Gen, dfm$ang_M), col="darkblue", lwd = 3) #add local regression
#   # lines(lowess(dfw$Gen, dfw$ang_M), col="yellowgreen", lwd = 3)
#   df500 <- subset(dfang1, Gen >= 500)
#   bymodel <- by(df500$ang_M, list(df500$Gen, df500$V10), FUN=mean.angle.pi)
#   lines(as.numeric(rownames(bymodel)), bymodel[,"m"], col="darkblue", lwd = 3)
#   lines(as.numeric(rownames(bymodel)), bymodel[,"w"], col="yellowgreen", lwd = 3)
#   lines(dfang1$Gen, dfang1$ang_S, col="orange", type = "l", lty=3, lwd = 2) #S orientation
#   axis(side=2, at=c(-pi, -pi/2, -pi/4, 0 , pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi), mgp = c(1.75, 0.75, 0))
#   legend <- FALSE
#   }
# dev.off()
# print("Ellipses done !")

#try to plot lines instead of points in genplot



