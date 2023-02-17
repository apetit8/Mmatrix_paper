source("scripts/functions_R/All_functions.R")
#########################################
sims.dirs <- list.dirs("simul/fig_2ab", recursive = FALSE)
angle = c(0.785, 0, -0.393)
of        <- "fig1"
modulo <- pi
#####################

M.factor <- 200 #400
G.factor <- 20
S.factor <- 1


#Df with generations
df.m.s.gen <- data.frame()
for (i in angle) {
  sims.dir  <- list.files(path=sims.dirs, pattern=paste0("simuangle", i,"$"), full.names=TRUE)
  df <- df.figsuppG(sims.dir, all.gen=TRUE)
  pop <- str_split(df$data.dir, "simul/fig_2ab/", n=2, simplify = TRUE)
  df[,10] <- sprintf("%s", pop[,2])
  model <- str_split(pop[,2], "/simuangle", n=2, simplify = TRUE)
  df[,11] <- sprintf("%s", model[,1])
  df.m.s.gen <- rbind(df.m.s.gen, df) 
}
df.m <- subset(df.m.s.gen, V11 == "m") #Df of Multilinear simulations
df.fkl <- subset(df.m.s.gen, V11 == "fkl")
df.w <- subset(df.m.s.gen, V11 == "w") #Df of Wagner simulations


colors <- c("maroon2", "darkblue", "yellowgreen")



# png(file="figures/fig1_part1.png", width=400, height=400)
cairo_pdf("figures/fig_supp3_a.pdf", width=6, height=6)
sims.dir  <- list.files(path=sims.dirs, pattern=paste0("simuangle", 0.785,"$"), full.names=TRUE)
# mar=c(0,0,0,0)
dfang1 <- subset(df.m.s.gen, round(ang_S, 1) == round(0.785, 1))
dfm <- subset(df.m, round(ang_S, 1) == round(0.785, 1))
dfw <- subset(df.w, round(ang_S, 1) == round(0.785, 1))
plot(dfang1$Gen, dfang1$ang_G, ylim =c(-pi/2,pi/2), xlim = c(min(dfang1$Gen), max(dfang1$Gen)),
     main="", yaxt="n", ylab = expression(paste("G matrix direction ",alpha, "(G)")), xlab = "Generation",
     col=alpha(colors[factor(dfang1$V11)],0.2),mar=c(1, 1, 0, 1), mgp = c(1.75, 0.75, 0),cex.lab=1.5, cex.axis=1.5)
df500 <- subset(dfang1, Gen >= 500)
bymodel <- by(df500$ang_G, list(df500$Gen, df500$V11), FUN=mean.angle.pi)
lines(as.numeric(rownames(bymodel)), bymodel[,"m"], col="darkblue", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"fkl"], col="maroon2", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"w"], col="yellowgreen", lwd = 3)
lines(dfang1$Gen, dfang1$ang_S, col="orange", type = "l", lty=3, lwd = 2) #S orientation
axis(side=2, at=c(-pi, -pi/2, -pi/4, 0 , pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi), mgp = c(1.75, 0.75, 0))
legend("bottomleft", lty=1, box.lty=0,  bg="transparent", col=c("yellowgreen","darkblue","darkred","orange"),
       legend=c(paste0("G GRN"), paste0("G multilinear"), paste0("G FKL"),paste0("S")))


oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, another_plot=TRUE,
                   xlim=c(-7,7), ylim=c(-7,7), all.reps=FALSE, xlab="", ylab="",xcoord=c(1000, 10000),ycoord=c(-2.5, 0.5),
                   yaxt = "n", xaxt = "n", legend=FALSE, asp=1, axes = FALSE, Gell=TRUE, Mell=FALSE,  mgp = c(0, 0, 0))
dev.off()
########

cairo_pdf("figures/fig_supp3_b.pdf", width=6, height=6)
sims.dir  <- list.files(path=sims.dirs, pattern=paste0("simuangle", -0.393,"$"), full.names=TRUE)
dfang1 <- subset(df.m.s.gen, round(ang_S, 1) == round(-0.393, 1))
dfm <- subset(df.m, round(ang_S, 1) == round(-0.393, 1))
dfw <- subset(df.w, round(ang_S, 1) == round(-0.393, 1))
plot(dfang1$Gen, dfang1$ang_G, ylim =c(-pi/2,pi/2), xlim = c(min(dfang1$Gen), max(dfang1$Gen)),
     main="", yaxt="n", ylab = expression(paste("G matrix direction ",alpha, "(G)")), xlab = "Generation", col=alpha(colors[factor(dfang1$V11)],0.2),
     mar=c(1, 1, 0, 1), mgp = c(1.75, 0.75, 0),cex.lab=1.5, cex.axis=1.5)
df500 <- subset(dfang1, Gen >= 500)
bymodel <- by(df500$ang_G, list(df500$Gen, df500$V11), FUN=mean.angle.pi)
lines(as.numeric(rownames(bymodel)), bymodel[,"m"], col="darkblue", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"fkl"], col="maroon2", lwd = 3)
lines(as.numeric(rownames(bymodel)), bymodel[,"w"], col="yellowgreen", lwd = 3)
lines(dfang1$Gen, dfang1$ang_S, col="orange", type = "l", lty=3, lwd = 2) #S orientation
axis(side=2, at=c(-pi, -pi/2, -pi/4, 0 , pi/4, pi/2, pi), labels=expression(-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi), mgp = c(1.75, 0.75, 0))
oneplot.allellipse(sims.dir, G.factor=G.factor, S.factor=S.factor, M.factor=M.factor, another_plot=TRUE,
                   xlim=c(-7,7), ylim=c(-7,7), all.reps=FALSE, xlab="", ylab="",xcoord=c(0, 10000),ycoord=c(-1, 2),
                   yaxt = "n", xaxt = "n", legend=FALSE, asp=1, axes = FALSE, Gell=TRUE, Mell=FALSE,  mgp = c(0, 0, 0))
dev.off()


################################################################################
sims.dirs <- list.dirs("simul/fig_2cd", recursive = FALSE)

  df.G <- data.frame(NULL)
  for (i in sims.dirs) {
    sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
    df <- df.figsuppG(sims.dir, all.gen=FALSE)
    df.G <- rbind(df.G, df)
  }
  df.G[,10] <- str_split(df.G$data.dir, "/", simplify = TRUE)[,3]


netw_names <- as_labeller(c(
  `0-grn` = "GRN model",
  `2-fkl` = "GP model",
  `1-mult` = "Multilinear Model"
))
#With eccentricity
pfig1de <- ggplot(data=df.G, aes(ang_S, ang_G))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  geom_point(aes(y=ang_G, col=ecc_G), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_G+pi, col=ecc_G), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_G-pi, col=ecc_G), alpha=0.2, show.legend = FALSE)+
  labs(y=expression(paste("Direction of mutational effects ",alpha, "(M)")), x=expression(paste("Fitness function direction, ",alpha, "(S)")), fill = expression("\u03BE\u03B1"))+
  scale_color_viridis_c(option = "plasma", limits=c(0,1), breaks = c(0.2, 0.5, 0.8))+
  labs(col = "M Eccentricity\ne(M)")+
  scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig1de <- pfig1de + facet_wrap(V10 ~., labeller = as_labeller(netw_names),  ncol=3) + theme(strip.background = element_blank())+ theme_bw(base_size = 13)+theme(panel.spacing = unit(0.7, "lines")) #, strip.text = element_blank()


netw_names <- as_labeller(c(
  `0-grn` = paste0("GRN model, \u03B2 = ", round(coef(lm(subset(df.G, V10=="0-grn")$corrG~ subset(df.G, V10=="0-grn")$corrS))[2] ,3)),
  `2-fkl` = paste0("GP model, \u03B2 = ", round(coef(lm(subset(df.G, V10=="2-fkl")$corrG~ subset(df.G, V10=="2-fkl")$corrS))[2] ,3)),
  `1-mult` = paste0("Multilinear model, \u03B2 = ", round(coef(lm(subset(df.G, V10=="1-mult")$corrG~ subset(df.G, V10=="1-mult")$corrS))[2] ,3))
))

pp <- ggplot(df.G, aes(x = corrS, y=corrG))+
  geom_point(aes(col=ecc_G), alpha=0.2)+theme(strip.background = element_blank())+ theme_bw(base_size = 13)+
  scale_color_viridis_c(option = "plasma", limits=c(0,1), breaks = c(0.2, 0.5, 0.8))+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, size=0.4)+
  labs(y=expression(paste("Mutational effect correlation r(M)")), x=expression(paste("Fitness function correlation, r(S)")), col = "M Eccentricity\ne(M)")+
  coord_fixed(ratio = 1)+
  facet_wrap(V10 ~., labeller = as_labeller(netw_names),  ncol=3)+
  theme(plot.margin = margin(t=4,0,0,0, "lines"),legend.direction="horizontal", legend.position = c(0.5, 1.27))+theme(panel.spacing = unit(0.7, "lines"))


cairo_pdf("figures/fig_supp3_c.pdf", width=8, height=4)
grid.arrange(
  pfig1de,
  ncol = 1,
  nrow = 1,
  clip = FALSE
)
dev.off()

cairo_pdf("figures/fig_supp3_d.pdf", width=8, height=4)
grid.arrange(
  pp,
  ncol = 1,
  nrow = 1,
  clip = FALSE
)
dev.off()

