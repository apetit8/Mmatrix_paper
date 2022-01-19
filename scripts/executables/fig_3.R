source("../functions_R/All_functions.R")
########################################################################################################################
sims.dirs <- list.dirs(
  "../../simul/fig_3/rnul", recursive = FALSE
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
  "../../simul/fig_3/rneg", recursive = FALSE
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
  "../../simul/fig_3/rpos2", recursive = FALSE
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

pdfname   <- print("../../figures/fig3-1.pdf")
cairo_pdf(pdfname, width=9, height=8)
grid.arrange(
  g1,g3,g4,
  ncol = 2,
  nrow = 2,
  widths = c(1,1),
  clip = FALSE
)
dev.off()

pdfname   <- print("../../figures/fig3.pdf")
cairo_pdf(pdfname, width=13, height=4.5)
grid.arrange(
  g1,g4,g3,
  ncol = 3,
  nrow = 1,
  widths = c(1,1,1.22),
  clip = FALSE
)
dev.off()

# library("DataCombine")
# dftest <- grepl.sub(data=df.all, pattern="FITNESS_OPTIMUM-0.5-FITNESS_OPTIMUM2-0.5", Var="V10", keep.found = TRUE)
# ggplot(data=dftest, aes(ang_S, ang_M, col=V11))+
#   coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
#   geom_abline()+
#   labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
#   geom_point()+
#   geom_smooth(alpha = 0, method = "AngleSmooth", span=0.9, level = 0.95)
# 
# mean(dftest$V11)
# mean(dftest$sq_dist)



######################
sims.dirs.1 <- list.dirs(
  "../../simul/fig_3/rround", recursive = FALSE
)
modulo <- pi
#####################

df.rround <- df.opt.map(sims.dirs.1)


g2 <- ggplot(df.rround, aes(x=P_mean_A, y=P_mean_B, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako", limits=c(0, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

g22 <- ggplot(df.rround, aes(x=P_mean_A, y=P_mean_B, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako", limits=c(-pi/2, pi/2))+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle mean for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Mean \u03B1")
g22


pdfname   <- print("../../figures/fig3-2.pdf")
cairo_pdf(pdfname, width=9, height=4)
grid.arrange(
  g2,g22,
  ncol = 2,
  nrow = 1,
  widths = c(1,1),
  clip = FALSE
)
dev.off()
