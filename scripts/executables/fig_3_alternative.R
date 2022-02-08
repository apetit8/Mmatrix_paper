source("../functions_R/All_functions.R")
###########################################
sims.dirs <- "../../simul/fig_3_bis/no_ab"
modulo <- pi
#####################

df.rnul <- df.opt.map2(sims.dirs, modulo=modulo)

g1 <- ggplot(df.rnul, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = FALSE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ S corr applied once phenotypes at steady state"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))
g1

############################################
sims.dirs <- list.dirs(
  "../../simul/fig_3/rnul", recursive = FALSE
)
modulo <- pi
#####################

df.rnul2 <- df.opt.map(sims.dirs, modulo=modulo)

g2 <- ggplot(df.rnul2, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.02, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("B/  Selectional correlation applied from the start"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))
g2

pdfname   <- print("../../figures/fig3_ter.pdf")
cairo_pdf(pdfname, width=10, height=6)
grid.arrange(
  g1, g2,
  ncol = 2,
  nrow = 1,
  widths = c(0.84,1),
  
  clip = FALSE
)
dev.off()

