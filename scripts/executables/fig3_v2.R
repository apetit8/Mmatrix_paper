source("../functions_R/All_functions.R")
#####################
sims.dirs <- c(
  "../../simul/fig_3_v2/no_ab"
)
modulo <- pi
#####################

dfv2.rnul <- df.opt.map2(sims.dirs, modulo=modulo, gen = FALSE)

dir <- str_split(dfv2.rnul$data.dir, ".par-R", n=3, simplify = TRUE)
dfv2.rnul$dir2 <- dir[,1]

dfv2_angevolv <- data.frame()
j <-1
for (u in unique(as.list(dir[,1]))){
  dfv2f <- subset(dfv2.rnul, dir2 == u)
  dfv2f[12] <- mean(dfv2f$sq_dist)
  dfv2_angevolv <- rbind(dfv2_angevolv, dfv2f[12] )
  j <-j+1
}
dfv2.rnul$d_alpha <- dfv2_angevolv[,1] #d_alpha

dir <- str_split(dfv2.rnul$dir2, ".FILE", n=3, simplify = TRUE)
dir <- str_split(dir[,1], "OPTIMUM-", n=3, simplify = TRUE)
dfv2.rnul$opti <- as.numeric(dir[,2]) #Optimum

# write.csv(dfv2.rnul, "../data/fig3-no_ab_fullgen_fullnetw.csv")

dfv2.rnul5000 <- subset(dfv2.rnul, Gen==5000)

dfv2.rnul <- dfv2.rnul5000
  
g1 <- ggplot(data=dfv2.rnul, aes(x=opti, y=ang_M) )+
  geom_point(aes(col=ang_S))+#, show.legend = FALSE
  stat_summary(fun=mean, colour="red", geom="line")+
  scale_color_gradient2(low="blue", mid="yellow",
                        high="blue", space ="Lab" )+
  labs(y="M direction", x="OPTIMUM")+
  theme_light()
g1

g2 <- ggplot(data=dfv2.rnul, aes(A_B, ang_M) )+
  geom_point(aes(col=B_A))+#, show.legend = FALSE
  # stat_summary(fun=mean, colour="red", geom="line")+
  scale_color_gradient2(low="blue", mid="yellow",
                        high="green", space ="Lab" )+
  labs(y="A_B", x="M direction")+
  theme_light()

g3 <- ggplot(data=dfv2.rnul, aes(A_B, B_A) )+
  geom_point(aes(col=opti))+#, show.legend = FALSE
  # stat_summary(fun=mean, colour="red", geom="line")+
  scale_color_gradient2(low="blue", mid="yellow",
                        high="green", space ="Lab" )+
  theme_light()

ggplot(data=dfv2.rnul, aes(A_B, B_A) )+
  geom_point(aes(col=ang_M))+#, show.legend = FALSE
  # stat_summary(fun=mean, colour="red", geom="line")+
  scale_color_gradient2(low="blue", mid="yellow",
                        high="green", space ="Lab" )+
  theme_light()


pfig3v2 <- ggplot(data=dfv2.rnul, aes(ang_S, ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  #previously : alpha of 0.1
  geom_point(aes(y=ang_M, col=B_A), alpha=0.16, show.legend = TRUE)+
  scale_colour_viridis_c(option = "plasma")+
  labs(y="M direction", x="S direction", fill = expression("\u03BE\u03B1"))
#scale_color_discrete(name=expression("\u03BE\u03B1"), labels=c("0.676", "0.403", "0.025","0.985", "0.505","0.373"))
pfig3v2 <- pfig3v2 + scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                                        labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig3v2 <- pfig3v2 + facet_wrap(vars(opti), ncol=3) + theme(strip.background = element_blank()) #, strip.text = element_blank()
pfig3v2

cairo_pdf("../../figures/fig_3_opti_and_things.pdf", width=17, height=5)
grid.arrange(
  g1, g2, g3,
  ncol = 3,
  nrow = 1,
  widths = c(1,1,1),
  clip = FALSE
)
dev.off()

################################################################################



















#################FIG for both



#To show : with some phenotypes, it is harder to have a good evolvability 
#BUT : not because there isn't any network for that, just because the evolutionary distance
#is too great. So actually, gene expression isn't enough to assume a correlation/ co-variance direction.


###########################################################

#Tab of correlation between the M direction and the W cell
tab_plast <- data.frame()
tab_plast[1:4,1:4] <- c(
  c(cor(dfv2.rnul$ang_M,dfv2.rnul$A_A), cor(dfv2.rnul$ang_M,dfv2.rnul$B_A), cor(dfv2.rnul$ang_M,dfv2.rnul$C_A), cor(dfv2.rnul$ang_M,dfv2.rnul$D_A)),
  c(cor(dfv2.rnul$ang_M,dfv2.rnul$A_B), cor(dfv2.rnul$ang_M,dfv2.rnul$B_B), cor(dfv2.rnul$ang_M,dfv2.rnul$C_B), cor(dfv2.rnul$ang_M,dfv2.rnul$D_B)),
  c(cor(dfv2.rnul$ang_M,dfv2.rnul$A_C), cor(dfv2.rnul$ang_M,dfv2.rnul$B_C), cor(dfv2.rnul$ang_M,dfv2.rnul$C_C), cor(dfv2.rnul$ang_M,dfv2.rnul$D_C)),
  c(cor(dfv2.rnul$ang_M,dfv2.rnul$A_D), cor(dfv2.rnul$ang_M,dfv2.rnul$B_D), cor(dfv2.rnul$ang_M,dfv2.rnul$C_D), cor(dfv2.rnul$ang_M,dfv2.rnul$D_D)) #Values from Jallet phd thesis
)
setnames(tab_plast, 1:4, c("A","B","C","D"))
tab_plast
g <- tableGrob(as.matrix(tab_plast),rows =  c("A","B","C","D"))
plot(g)

#Tab of correlation between the optimum and the W cell
tab_plast <- data.frame()
tab_plast[1:4,1:4] <- c(
  c(cor(dfv2.rnul$opti,dfv2.rnul$A_A), cor(dfv2.rnul$opti,dfv2.rnul$B_A), cor(dfv2.rnul$opti,dfv2.rnul$C_A), cor(dfv2.rnul$opti,dfv2.rnul$D_A)),
  c(cor(dfv2.rnul$opti,dfv2.rnul$A_B), cor(dfv2.rnul$opti,dfv2.rnul$B_B), cor(dfv2.rnul$opti,dfv2.rnul$C_B), cor(dfv2.rnul$opti,dfv2.rnul$D_B)),
  c(cor(dfv2.rnul$opti,dfv2.rnul$A_C), cor(dfv2.rnul$opti,dfv2.rnul$B_C), cor(dfv2.rnul$opti,dfv2.rnul$C_C), cor(dfv2.rnul$opti,dfv2.rnul$D_C)),
  c(cor(dfv2.rnul$opti,dfv2.rnul$A_D), cor(dfv2.rnul$opti,dfv2.rnul$B_D), cor(dfv2.rnul$opti,dfv2.rnul$C_D), cor(dfv2.rnul$opti,dfv2.rnul$D_D)) #Values from Jallet phd thesis
)
setnames(tab_plast, 1:4, c("A","B","C","D"))
tab_plast
gg <- tableGrob(as.matrix(tab_plast),rows =  c("A","B","C","D"))
plot(gg)



#Tab of correlation between the s direction and the W cell
tab_plast <- data.frame()
tab_plast[1:4,1:4] <- c(
  c(cor(dfv2.rnul$ang_S,dfv2.rnul$A_A), cor(dfv2.rnul$ang_S,dfv2.rnul$B_A), cor(dfv2.rnul$ang_S,dfv2.rnul$C_A), cor(dfv2.rnul$ang_S,dfv2.rnul$D_A)),
  c(cor(dfv2.rnul$ang_S,dfv2.rnul$A_B), cor(dfv2.rnul$ang_S,dfv2.rnul$B_B), cor(dfv2.rnul$ang_S,dfv2.rnul$C_B), cor(dfv2.rnul$ang_S,dfv2.rnul$D_B)),
  c(cor(dfv2.rnul$ang_S,dfv2.rnul$A_C), cor(dfv2.rnul$ang_S,dfv2.rnul$B_C), cor(dfv2.rnul$ang_S,dfv2.rnul$C_C), cor(dfv2.rnul$ang_S,dfv2.rnul$D_C)),
  c(cor(dfv2.rnul$ang_S,dfv2.rnul$A_D), cor(dfv2.rnul$ang_S,dfv2.rnul$B_D), cor(dfv2.rnul$ang_S,dfv2.rnul$C_D), cor(dfv2.rnul$ang_S,dfv2.rnul$D_D)) #Values from Jallet phd thesis
)
setnames(tab_plast, 1:4, c("A","B","C","D"))
tab_plast
ggg <- tableGrob(as.matrix(tab_plast),rows =  c("A","B","C","D"))
plot(ggg)



#Tab of correlation between the s direction and the W cell
tab_plast <- data.frame()
tab_plast[1:4,1:4] <- c(
  c(cor(dfv2.rnul5000$opti,dfv2.rnul5000$A_A), cor(dfv2.rnul5000$opti,dfv2.rnul5000$B_A), cor(dfv2.rnul5000$opti,dfv2.rnul5000$C_A), cor(dfv2.rnul5000$opti,dfv2.rnul5000$D_A)),
  c(cor(dfv2.rnul5000$opti,dfv2.rnul5000$A_B), cor(dfv2.rnul5000$opti,dfv2.rnul5000$B_B), cor(dfv2.rnul5000$opti,dfv2.rnul5000$C_B), cor(dfv2.rnul5000$opti,dfv2.rnul5000$D_B)),
  c(cor(dfv2.rnul5000$opti,dfv2.rnul5000$A_C), cor(dfv2.rnul5000$opti,dfv2.rnul5000$B_C), cor(dfv2.rnul5000$opti,dfv2.rnul5000$C_C), cor(dfv2.rnul5000$opti,dfv2.rnul5000$D_C)),
  c(cor(dfv2.rnul5000$opti,dfv2.rnul5000$A_D), cor(dfv2.rnul5000$opti,dfv2.rnul5000$B_D), cor(dfv2.rnul5000$opti,dfv2.rnul5000$C_D), cor(dfv2.rnul5000$opti,dfv2.rnul5000$D_D)) #Values from Jallet phd thesis
)
setnames(tab_plast, 1:4, c("A","B","C","D"))
tab_plast
ggg <- tableGrob(as.matrix(tab_plast),rows =  c("A","B","C","D"))
plot(ggg)


