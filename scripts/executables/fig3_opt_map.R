source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/ggplot2_custom_functions.R")
########################################################################################################################
sims.dirs <- c(
  "../../simul/fig_3/raster_map3/0pi", "../../simul/fig_3/raster_map2/0pi", "../../simul/fig_3/raster_map2/-pi_4",  "../../simul/fig_3/raster_map2/pi_4",
  "../../simul/fig_3/raster_map2/-3pi_8","../../simul/fig_3/raster_map2/3pi_8","../../simul/fig_3/raster_map2/pi_8",
  "../../simul/fig_3/raster_map2/-pi_8",  "../../simul/fig_3/raster_map3/-pi_4" ,  "../../simul/fig_3/raster_map3/pi_4","../../simul/fig_3/raster_map3/-3pi_8","../../simul/fig_3/raster_map3/3pi_8",
  "../../simul/fig_3/raster_map3/pi_8","../../simul/fig_3/raster_map3/-pi_8","../../simul/fig_3/raster_map3/1.57.temp"
)
of        <- "fig3"
modulo <- pi
#####################

df.all <- data.frame(NULL)
for (i in sims.dirs) {
  df <- df.raster.map(i, modulo=modulo)
  pop <- str_split(i, "../../simul/fig_3/raster_map", n=2, simplify = TRUE)
  df[,9] <- sprintf("%s", pop[,2])
  df.all <- rbind(df.all, df)
}

g1 <- ggplot(df.all, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = "mean")+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.013, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("B/ Alignment score for different optimum\nphenotypes with correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression(Delta*a))


######################
sims.dirs.1 <- list.dirs(
  "../../simul/fig_3/round_map", recursive = FALSE
)
modulo <- pi
#####################

df.all1 <- data.frame(NULL)
for (i in sims.dirs.1) {
  df <- df.raster.map(i, modulo=modulo)
  pop <- str_split(i, "../../simul/fig_3/round_map", n=2, simplify = TRUE)
  df[,9] <- sprintf("%s", pop[,2])
  df.all1 <- rbind(df.all1, df)
}

# g2 <- ggplot(df.all1, aes(x=P_mean_A, y=P_mean_B, z=ang_M))+
#   coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
#   stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
#   scale_fill_viridis_c(option = "magma",limits=c(-pi/2,pi/2))+
#   geom_point(aes(y=0, x=0))+
#   labs(title = ("B/ Mean M angle for different optimum phenotypes"), x="Expression gene A", y="Expression gene B")+
#   labs(fill = "M angle")

# ggplot(df.all1, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
#   coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
#   stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = "mean")+
#   scale_fill_viridis_c(option = "plasma",limits=c(-0.08, 0.89))+
#   geom_point(aes(y=0, x=0))+
#   labs(title = ("B/ Mean delta a for different optimum phenotypes"), x="Expression gene A", y="Expression gene B")+
#   labs(fill = "delta a")

g2 <- ggplot(df.all1, aes(x=P_mean_A, y=P_mean_B, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = st_dev_abs)+
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

df.all5 <- data.frame(NULL)
for (i in sims.dirs.5) {
  df <- df.raster.map(i, modulo=modulo)
  pop <- str_split(i, "../../simul/fig_3/raster_initcorrneg", n=2, simplify = TRUE)
  df[,9] <- sprintf("%s", pop[,2])
  df.all5 <- rbind(df.all5, df)
}

g3 <-   ggplot(df.all5, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = "mean")+
  # scale_fill_viridis_c(option = "plasma")+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.013, 1))+
  labs(title =  ("C/ Alignment score for different optimum\nphenotypes with a negative initial correlation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression(Delta*a))

######################
sims.dirs.6 <- list.dirs(
  "../../simul/fig_3/raster_initcorrpos", recursive = FALSE
)
modulo <- pi
#####################

df.all6 <- data.frame(NULL)
for (i in sims.dirs.6) {
  df <- df.raster.map(i, modulo=modulo)
  pop <- str_split(i, "../../simul/fig_3/raster_initcorrpos", n=2, simplify = TRUE)
  df[,9] <- sprintf("%s", pop[,2])
  df.all6 <- rbind(df.all6, df)
}

g4 <- ggplot(df.all6, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = "mean")+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.013, 1))+
  labs(title =  ("D/ Alignment score for different optimum\nphenotypes with a positive initial correlation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression(Delta*a))


#################################
# PDF
#################################

pdfname   <- print(sprintf("../../figures/%s_selection_param.pdf", of))
pdf(pdfname, width=9, height=8)

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

