source("scripts/functions_R/All_functions.R")
library(png)
#####################
sims.dirs <- list.dirs("simul/fig_2cd", recursive = FALSE)
#####################

#Data
df.fig2cd <- df.data(sims.dirs, pattern = "simul/fig_2cd/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE)

netw_names <- as_labeller(c(
  `0-grn` = "GRN model",
  `2-fkl` = "GP model",
  `1-mult` = "Multilinear Model"
))

#With eccentricity
pfig2cd <- ggplot(data=df.fig2cd, aes(ang_S, ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  geom_point(aes(y=ang_M, col=ecc_M), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, col=ecc_M), alpha=0.2, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, col=ecc_M), alpha=0.2, show.legend = FALSE)+
  labs(y=expression(paste("Direction of mutational effects ",alpha, "(M)")), x=expression(paste("Fitness function direction, ",alpha, "(S)")), fill = expression("\u03BE\u03B1"))+
  scale_color_viridis_c(option = "plasma")+
  labs(col = "M Eccentricity\ne(M)")+
  scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig2cd <- pfig2cd + facet_wrap(pop ~., labeller = as_labeller(netw_names),  ncol=3) + theme(strip.background = element_blank())+ theme_bw(base_size = 13) #base_size = 12


cairo_pdf("figures/fig_2c.pdf", width=8, height=4)
grid.arrange(
  pfig2cd,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
dev.off()


netw_names <- as_labeller(c(
  `0-grn` = paste0("GRN model, \u03B2 = ", round(coef(lm(subset(df.fig2cd, pop=="3-grn")$corrM~ subset(df.fig2cd, pop=="3-grn")$corrS))[2] ,3)),
  `2-fkl` = paste0("GP model, \u03B2 = ", round(coef(lm(subset(df.fig2cd, pop=="2-fkl")$corrM~ subset(df.fig2cd, pop=="2-fkl")$corrS))[2] ,3)),
  `1-mult` = paste0("Multilinear model, \u03B2 = ", round(coef(lm(subset(df.fig2cd, pop=="1-mult")$corrM~ subset(df.fig2cd, pop=="1-mult")$corrS))[2] ,3))
))

pp <- ggplot(df.fig2cd, aes(x = corrS, y=corrM))+
  geom_point(aes(col=ecc_M), alpha=0.2)+theme(strip.background = element_blank())+theme_bw(base_size = 13)+
  scale_color_viridis_c(option = "plasma")+
  labs(y=expression(paste("Mutational effect correlation r(M)")), x=expression(paste("Fitness function correlation, r(S)")), col = "M Eccentricity\ne(M)")+
  coord_fixed(ratio = 1)+
  facet_wrap(pop ~., labeller = as_labeller(netw_names),  ncol=3)+
  theme(plot.margin = margin(t=4,0,0,0, "lines"),legend.direction="horizontal", legend.position = c(0.5, 1.27))

cairo_pdf("figures/fig_2d.pdf", width=8, height=4)
grid.arrange(
  pp,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
dev.off()



