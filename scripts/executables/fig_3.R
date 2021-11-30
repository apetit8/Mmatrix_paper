source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/ggplot2_custom_functions.R")
########################################################################################################################
sims.dirs <- list.dirs(
  "../../simul/fig_3/raster_initcorrnul", recursive = FALSE
)
of        <- "fig3"
modulo <- pi
#####################

df.all <- df.opt.map(sims.dirs, modulo=modulo)

g1 <- ggplot(df.all, aes(x=P_mean_A, y=P_mean_B, z=V11))+
  coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = "mean")+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.05, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("B/ Alignment score for different optimum\nphenotypes with correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))

######################
sims.dirs.1 <- list.dirs(
  "../../simul/fig_3/round_map", recursive = FALSE
)
modulo <- pi
#####################

df.all1 <- df.opt.map(sims.dirs.1)


g2 <- ggplot(df.all1, aes(x=P_mean_A, y=P_mean_B, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
  sstat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako")+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

######################
sims.dirs.5 <- list.dirs(
  "../../simul/fig_3/raster_initcorrneg", recursive = FALSE
)
modulo <- pi
#####################

df.all5 <- df.opt.map(sims.dirs.5)

g3 <-   ggplot(df.all5, aes(x=P_mean_A, y=P_mean_B, z=V11))+
  coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = "mean")+
  # scale_fill_viridis_c(option = "plasma")+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.05, 1))+
  labs(title =  ("C/ Alignment score for different optimum\nphenotypes with a negative initial correlation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))

######################
sims.dirs.6 <- list.dirs(
  "../../simul/fig_3/raster_initcorrpos", recursive = FALSE
)
modulo <- pi
#####################

df.all6 <- df.opt.map(sims.dirs.6)

g4 <- ggplot(df.all6, aes(x=P_mean_A, y=P_mean_B, z=V11))+
  coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = "mean")+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.05, 1))+
  labs(title =  ("D/ Alignment score for different optimum\nphenotypes with a positive initial correlation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))


#################################
# PDF
#################################

pdfname   <- print(sprintf("../../figures/%s.pdf", of))
cairo_pdf(pdfname, width=9, height=8)

lay <- rbind(c(1,2),c(3,4))

grid.arrange(
  g2,g1,g3,g4,
  ncol = 2,
  nrow = 2,
  widths = c(1,1),
  layout_matrix = lay,
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
