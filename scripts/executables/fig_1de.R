source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs <- list.dirs("../../simul/fig_1de", recursive = FALSE)
modulo <- pi
#####################

#Data
df.fig1de <- df.data(sims.dirs, pattern = "../../simul/fig_1de/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE)

netw_names <- as_labeller(c(
  `2-full` = "GRN model",
  `1-mult` = "Multilinear Model"
))

#With eccentricity
pfig1de <- ggplot(data=df.fig1de, aes(ang_S, ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  geom_point(aes(y=ang_M, col=ecc_M), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, col=ecc_M), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, col=ecc_M), alpha=0.2, show.legend = TRUE)+
  labs(y=expression(paste("Direction of mutational effects ",alpha, "(M)")), x=expression(paste("Fitness function, ",alpha, "(S)")), fill = expression("\u03BE\u03B1"))+
  scale_color_viridis_c(option = "plasma")+
  labs(col = "M Eccentricity e(M)")+
  scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig1de <- pfig1de + facet_wrap(pop ~., labeller = as_labeller(netw_names),  ncol=3) + theme(strip.background = element_blank())+ theme_bw() #, strip.text = element_blank()
pfig1de

cairo_pdf("../../figures/fig1_de.pdf", width=7, height=4)
grid.arrange(
  pfig1de,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
dev.off()
