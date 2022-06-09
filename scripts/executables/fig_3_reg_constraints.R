source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs <- list.dirs("../../simul/fig_3", recursive = FALSE)
modulo <- pi
#####################

#Data
df.fig2 <- df.data(sims.dirs, pattern = "../../simul/fig_3/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE)

netw_names <- as_labeller(c(
  `3-full` = "C/ Free full network",
  `4-full_-0.5` = "D/ AB regulation fixed at -0.5",
  `5-full_0` = "E/ AB regulation fixed at 0",
  `6-full_0.5` = "F/ AB regulation fixed at 0.5",
  `1-mult` = "A/ Multilinear Model",
  `2-sub_ntw_ac_bd` = "B/ A and B in different subnetwork"
))

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
pfig2 <- pfig2 + facet_wrap(pop ~., labeller = as_labeller(netw_names),  ncol=3) + theme(strip.background = element_blank())+ theme_bw() #, strip.text = element_blank()
pfig2

cairo_pdf("../../figures/fig_3_reg_constraints.pdf", width=10, height=7)
grid.arrange(
  pfig2,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
dev.off()




pfig2 <- ggplot(data=df.fig2, aes(ang_S, ang_M,col=netw))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  #previously : alpha of 0.1
  geom_point(aes(y=ang_M, col=pop), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, col=pop), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, col=pop), alpha=0.16, show.legend = TRUE)+
  # geom_point(aes(y=mean_ang_M, fill=netw), show.legend = FALSE)+
  labs(y="M direction", x="S direction", fill = expression("\u03BE\u03B1"))+
  scale_color_discrete(name=expression("\u03BE\u03B1"),
                       labels=c("0.676", "0.403", "0.025","0.985", "0.505","0.373"))

pfig2 <- pfig2 + scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                                    labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig2 <- pfig2 + facet_wrap(pop ~., labeller = as_labeller(netw_names),  ncol=3) + theme(strip.background = element_blank())+ theme_bw() #, strip.text = element_blank()
pfig2


cairo_pdf("../../figures/fig_3_color.pdf", width=10, height=7)
grid.arrange(
  pfig2,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
dev.off()







#####################
sims.dirs <- list.dirs("../../simul/witness", recursive = FALSE)
modulo <- pi
#####################

#Data
df.fig2 <- df.data(sims.dirs, pattern = "../../simul/witness/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE)

pfig2 <- ggplot(data=df.fig2, aes(ang_S, ang_M,col=netw))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  #previously : alpha of 0.1
  geom_point(aes(y=ang_M, col=pop), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, col=pop), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, col=pop), alpha=0.16, show.legend = TRUE)+
  # geom_point(aes(y=mean_ang_M, fill=netw), show.legend = FALSE)+
  labs(y="M direction", x="S direction", fill = expression("\u03BE\u03B1"))+
  scale_color_discrete(name=expression("\u03BE\u03B1"),
                       labels=c("0.676", "0.403", "0.025","0.985", "0.505","0.373"))

pfig2 <- pfig2 + scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                                    labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig2 <- pfig2 + facet_wrap(pop ~., labeller = as_labeller(netw_names),  ncol=3) + theme(strip.background = element_blank())+ theme_bw() #, strip.text = element_blank()
pfig2

cairo_pdf("../../figures/fig_witness.pdf", width=4, height=5)
grid.arrange(
  pfig2,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
dev.off()