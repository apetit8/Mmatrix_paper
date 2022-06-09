source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs.fig2c <-list.dirs(
  "../../simul/fig_4", recursive = FALSE
)
modulo <- pi
#####################

#Data
df.fig2c <- df.data(sims.dirs.fig2c, pattern = "../../simul/fig_4/", variable="netw")

# dir <- str_split(df.fig2c$data.dir, "simul", n=3, simplify = TRUE)
# df.fig2c[,11] <- dir[,2]
# df_angevolv <- data.frame()
# mean_diff <- data.frame()
# j <-1
# for (u in as.list(df.fig2c[,11])){
#   dist <-  1- mse.ms(subset(df.fig2c, V11 == u)$ang_diff)/((pi^2)/12)
#   df_angevolv[j,1] <- dist
#   mean_diff[j,1] <- mean(subset(df.fig2c, V11 == u)$ang_diff)
#   j <-j+1
# }
# df.fig2c$d_alpha <- df_angevolv[,1]
# df.fig2c$m_diff <- mean_diff[,1]
# df.fig2c$dist2 <- 1- df.fig2c$ang_diff/((pi^2)/12)

ggplot(data=df.fig2c, aes(ang_S, d_alpha) )+
  geom_point(aes(col=netw ))+
  geom_line(aes(col=netw ))+
  # geom_smooth(aes(col=netw ), alpha=0, span=0.2)+
  coord_polar(start = pi/2, direction = -1)+
  scale_x_continuous(limits=c(-pi,pi),breaks=seq(-pi/2, pi/2, pi/4),labels=c("-\u03c0/2", "-\u03c0/4", "0", "\u03c0/4", "\u03c0/2"))+
  labs(x="S direction", y="d\u03B1")+
  scale_color_discrete(name = "A and B regulation", labels = c("Direct and Indirect", "Direct only", "Indirect : 1 gene", "Indirect: 2 genes", "Indirect: 3 genes"))+
theme_bw()#+  scale_color_manual(values=c("#F8766D", "#00BFC4", "#7CAE00","#C77CFF"))



#################Circular plot ?

  
  pfig2 <- ggplot(data=df.fig2c, aes(ang_S, ang_M))+
    coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
    geom_abline(colour="#666666")+
    geom_abline(intercept=pi, colour="#666666")+
    geom_abline(intercept=-pi, colour="#666666")+
    #previously : alpha of 0.1
    geom_point(aes(y=ang_M, col=1-ecc_M), alpha=0.16, show.legend = FALSE)+
    geom_point(aes(y=ang_M_mpi, col=1-ecc_M), alpha=0.16, show.legend = FALSE)+
    geom_point(aes(y=ang_M_ppi, col=1-ecc_M), alpha=0.16, show.legend = TRUE)+
    # geom_point(aes(y=mean_ang_M, fill=netw), show.legend = FALSE)+
    labs(y="M direction", x="S direction", fill = expression("\u03BE\u03B1"))+
    scale_color_viridis_c(option = "plasma")+
    labs(col = "M Eccentricity")+
    scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                                      labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
    scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                       labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
  pfig2 <- pfig2 + facet_wrap(vars(pop), ncol=3) + theme(strip.background = element_blank())+ theme_bw()  #, strip.text = element_blank()
  pfig2
  
g_sdfig2 <- table.index(sims.dirs.fig2, pattern = "../../simul/fig_2_bis/", ref = df.fig2c, bymean = FALSE, asgrob=TRUE)
plot(g_sdfig2)


##

df.fig2c$moduloquart <- modulo.all(df.fig2c$ang_M, modulo=pi/2)

lastplot <- ggplot(df.fig2c, aes(x=netw, y=moduloquart, col=d_alpha))+
  geom_boxplot()+
  geom_point()+ #shape=1
  scale_color_viridis_c(option = "plasma")+
  scale_y_continuous(limits=c(-pi/4,pi/4),breaks=seq(-pi/4, pi/4, pi/8),
                     labels=c("-\u03c0/4", "-\u03c0/8", "0", "\u03c0/8", "\u03c0/4"))+
  scale_x_discrete(labels = c("Direct and Indirect", "Direct only", "Indirect : 1 gene", "Indirect: 2 genes", "Indirect: 3 genes"))+
  labs(x="Regulation between selected genes", y="M direction", col=expression("d\u03B1"))+ theme_bw()


cairo_pdf("../../figures/fig_2_distrib.pdf", width=8, height=5)
grid.arrange(
  lastplot,
  ncol = 1,
  nrow = 1,
  clip = FALSE
)
dev.off()



#Ver 1
pfig2 <- ggplot(data=df.fig2c, aes(ang_S, ang_M,col=pop))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  #previously : alpha of 0.1
  geom_point(aes(y=ang_M, col=pop), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, col=pop), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, col=pop), alpha=0.16, show.legend = TRUE)+
  # geom_point(aes(y=mean_ang_M, fill=netw), show.legend = FALSE)+
  labs(y="M direction", x="S direction", fill = expression("\u03BE\u03B1"))
pfig2 <- pfig2 + scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                                    labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig2 <- pfig2 + facet_wrap(vars(pop), ncol=3) + theme(strip.background = element_blank()) #, strip.text = element_blank()
pfig2

cairo_pdf("../../figures/fig_2_ter3.pdf", width=9, height=9)
pfig2
dev.off()
