source("../functions_R/All_functions.R")
########################################################################################################################
sims.dirs <- list.dirs(
  "../../simul/fig_3_bis/no_ab", recursive = FALSE
)
of        <- "fig3"
modulo <- pi
#####################

df.rnul <- df.opt.map(sims.dirs, modulo=modulo)

g1 <- ggplot(df.rnul, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = FALSE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.02, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Without initial AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))

######################
sims.dirs.5 <- list.dirs(
  "../../simul/fig_3_bis/neg_ab", recursive = FALSE
)
modulo <- pi
#####################

df.rneg <- df.opt.map(sims.dirs.5)

g3 <-   ggplot(df.rneg, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean")+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.02, 1))+
  labs(title =  ("C/ Negative initial AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))
g3

######################
sims.dirs.6 <- list.dirs(
  "../../simul/fig_3_bis/pos_ab", recursive = FALSE
)
modulo <- pi
#####################

df.rpos<- df.opt.map(sims.dirs.6)

g4 <- ggplot(df.rpos, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = FALSE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.02, 1))+
  labs(title =  ("B/ Positive initial AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))
g4

#################################
# PDF
#################################

pdfname   <- print("../../figures/fig3bis.pdf")
cairo_pdf(pdfname, width=13, height=4.5)
grid.arrange(
  g1,g4,g3,
  ncol = 3,
  nrow = 1,
  widths = c(1,1,1.22),
  clip = FALSE
)
dev.off()
