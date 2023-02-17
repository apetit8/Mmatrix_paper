#Figure initW
source("scripts/functions_R/All_functions.R")
#####################
sims.dirs <- list.dirs("simul/fig_supp6_initW", recursive = FALSE)
#####################

#Data
df.initW <- df.data(sims.dirs, pattern = "simul/fig_supp6_initW/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE, all.gen=TRUE)
#Effet des mutations dans le réseau : pas toujours une lois normale ? (car pas la même moyenne retrouvée entre plusieurs tirages aléatoires)

netw_names <- as_labeller(c(
  `0` = "Generation 0",
  `10000` = "Generation 10000"
))

#With eccentricity
piw <- ggplot(data=subset(df.initW, Gen==0 | Gen==10000) , aes(ang_S, ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  geom_point(aes(y=ang_M, col=ecc_M), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, col=ecc_M), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, col=ecc_M), alpha=0.2, show.legend = TRUE)+
  labs(y=expression(paste(alpha, "(M)")), x=expression(paste(alpha, "(S)")))+
  scale_color_viridis_c(option = "plasma", limits=c(0.01,0.99))+
  labs(col = "e(M)")+
  scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
piw <- piw + facet_wrap(Gen ~., labeller = as_labeller(netw_names),  ncol=2) +
  theme_bw() #, strip.text = element_blank()

cairo_pdf("figures/fig_supp6_initW.pdf", width=8, height=4)
grid.arrange(
  piw,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
dev.off()

piw2 <- ggplot(subset(df.initW, Gen==0 | Gen==10000) , aes(x = corrS, y=corrM))+
  geom_point(aes(col=ecc_M), alpha=0.2)+theme(strip.background = element_blank())+theme_bw(base_size = 13)+
  scale_color_viridis_c(option = "plasma")+
  labs(y=expression(paste("Mutational effect correlation r(M)")), x=expression(paste("Fitness function correlation, r(S)")), col = "M Eccentricity\ne(M)")+
  coord_fixed(ratio = 1)+
  facet_wrap(Gen ~., labeller = as_labeller(netw_names),  ncol=3)+
  theme(plot.margin = margin(t=4,0,0,0, "lines"),legend.direction="horizontal", legend.position = c(0.5, 1.27))

cairo_pdf("figures/fig_supp6_initW_corr.pdf", width=8, height=4)
grid.arrange(
  piw2,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
dev.off()
