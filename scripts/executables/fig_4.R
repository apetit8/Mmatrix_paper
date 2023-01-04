source("scripts/functions_R/All_functions.R")
library(png)
#####################
sims.dirs <- list.dirs("simul/fig_4", recursive = FALSE)
#####################

#Data
df.fig4 <- df.data(sims.dirs, pattern = "simul/fig_4/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE)

netw_names <- as_labeller(c(
  `1-full_-0.5_fixed` = "a-b regulations fixed at -0.5",
  `3-full_0_fixed` = "a-b regulations fixed at 0",
  `2-full_0.5_fixed` = "a-b regulations fixed at 0.5"
))

#With eccentricity
pfig4 <- ggplot(data=df.fig3, aes(ang_S, ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  geom_point(aes(y=ang_M, col=ecc_M), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, col=ecc_M), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, col=ecc_M), alpha=0.2, show.legend = TRUE)+
  labs(y=expression(paste(alpha, "(M)")), x=expression(paste(alpha, "(S)")))+
  scale_color_viridis_c(option = "plasma")+
  labs(col = "e(M)")+
  scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig4 <- pfig4 + facet_wrap(pop ~., labeller = as_labeller(netw_names),  ncol=3) +
  theme_bw()+ theme(plot.margin = unit(c(3.3, 0, 0, 0), "cm")) #, strip.text = element_blank()



png.netw1 = readPNG('templates/fig_4/networks.png')

cairo_pdf("figures/fig_4.pdf", width=8, height=5.8)
grid.arrange(
  pfig4,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
grid.raster(png.netw1, x=0.49, y=0.76, width=0.68)
dev.off()

