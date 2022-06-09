source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs <- list.dirs("../../simul/fig_light_constraint2", recursive = FALSE)
modulo <- pi
#####################

#Data
df.fig2 <- df.data(sims.dirs, pattern = "../../simul/fig_light_constraint2/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE)

#With eccentricity
pfig2 <- ggplot(data=df.fig2, aes(ang_S, ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  geom_point(aes(y=ang_M, col=1-ecc_M), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, col=1-ecc_M), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, col=1-ecc_M), alpha=0.16, show.legend = TRUE)+
  labs(y="M direction", x="S direction", fill = expression("\u03BE\u03B1"))+
  scale_color_viridis_c(option = "plasma")+
  labs(col = "M Eccentricity")+
  scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig2 <- pfig2 + facet_wrap(pop ~., ncol=3) + theme(strip.background = element_blank())+ theme_bw() #, strip.text = element_blank()
pfig2