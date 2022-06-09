source("../functions_R/All_functions.R")
################################################################################
sims.dirs <- "../../simul/fig_3/no_ab"
of        <- "fig3"
modulo <- pi
#####################

df.rnul <- df.opt.map2(sims.dirs, modulo=modulo, gen=TRUE)
write.csv(df.rnul, "../data/fig3-no_ab_fullgen.csv")

dir <- str_split(df.rnul$data.dir, ".par-R", n=3, simplify = TRUE)
df.rnul[,11] <- dir[,1]
dir <- str_split(df.rnul$V11, ".FILE", n=3, simplify = TRUE)
dir <- str_split(dir[,1], "FITNESS_OPTIMUM-", n=3, simplify = TRUE)
dir <- str_split(dir[,2], "-FITNESS_OPTIMUM2-", n=3, simplify = TRUE)
df.rnul[,12] <- dir[,1]
df.rnul[,13] <- dir[,2]

df.rnul[,12] <- as.numeric(df.rnul[,12])
df.rnul[,13] <- as.numeric(df.rnul[,13])


df.rnul10000 <- subset(df.rnul, Gen==10000)
df.rnul5000 <- subset(df.rnul, Gen==5000)
df.rnul5250 <- subset(df.rnul, Gen==5250)

ggplot(df.rnul5250, aes(x=V12, y=V13, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako")+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

ggplot(df.rnul5000, aes(x=V12, y=V13, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako")+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

ggplot(df.rnul10000, aes(x=V12, y=V13, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako", limits=c(0, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

ggplot(df.rnul5000, aes(x=V12, y=V13, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))

ggplot(df.rnul10000, aes(x=V12, y=V13, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))


ggplot(df.rnul5000, aes(x=V12, y=V13, z=Fitness))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("Fitness"))

ggplot(df.rnul10000, aes(x=V12, y=V13, z=Fitness))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("Fitness"))

ggplot(df.rnul5250, aes(x=V12, y=V13, z=Fitness))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("Fitness"))

ggplot(df.rnul5250, aes(x=V12, y=V13, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))


######################
sims.dirs.5 <- "../../simul/fig_3/neg_ab"
modulo <- pi
#####################

df.rneg.allgen <- df.opt.map2(sims.dirs.5, gen=TRUE)
write.csv(df.rneg.allgen, "../data/fig3-neg_ab_fullgen.csv")

dir <- str_split(df.rneg.allgen$data.dir, ".par-R", n=3, simplify = TRUE)
df.rneg.allgen[,11] <- dir[,1]
dir <- str_split(df.rneg.allgen$V11, ".FILE", n=3, simplify = TRUE)
dir <- str_split(dir[,1], "FITNESS_OPTIMUM-", n=3, simplify = TRUE)
dir <- str_split(dir[,2], "-FITNESS_OPTIMUM2-", n=3, simplify = TRUE)
df.rneg.allgen[,12] <- dir[,1]
df.rneg.allgen[,13] <- dir[,2]

df.rneg.allgen[,12] <- as.numeric(df.rneg.allgen[,12])
df.rneg.allgen[,13] <- as.numeric(df.rneg.allgen[,13])

df.rneg10000 <- subset(df.rneg.allgen, Gen==10000)
df.rneg5000 <- subset(df.rneg.allgen, Gen==5000)
df.rneg5250 <- subset(df.rneg.allgen, Gen==5250)


ggplot(df.rneg5250, aes(x=V12, y=V13, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako")+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

ggplot(df.rneg5000, aes(x=V12, y=V13, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako")+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

ggplot(df.rneg10000, aes(x=V12, y=V13, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako", limits=c(0, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

ggplot(df.rneg5000, aes(x=V12, y=V13, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))

ggplot(df.rneg10000, aes(x=V12, y=V13, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))


ggplot(df.rneg5000, aes(x=V12, y=V13, z=Fitness))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("Fitness"))

ggplot(df.rneg10000, aes(x=V12, y=V13, z=Fitness))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("Fitness"))

ggplot(df.rneg5250, aes(x=V12, y=V13, z=Fitness))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("Fitness"))

ggplot(df.rneg5250, aes(x=V12, y=V13, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))



######################
sims.dirs.6 <- "../../simul/fig_3/pos_ab"
modulo <- pi
#####################

df.rpos.allgen<- df.opt.map2(sims.dirs.6, gen = TRUE)
write.csv(df.rpos.allgen, "../data/fig3-pos_ab_fullgen.csv")


dir <- str_split(df.rpos.allgen$data.dir, ".par-R", n=3, simplify = TRUE)
df.rpos.allgen[,11] <- dir[,1]
dir <- str_split(df.rpos.allgen$V11, ".FILE", n=3, simplify = TRUE)
dir <- str_split(dir[,1], "FITNESS_OPTIMUM-", n=3, simplify = TRUE)
dir <- str_split(dir[,2], "-FITNESS_OPTIMUM2-", n=3, simplify = TRUE)
df.rpos.allgen[,12] <- dir[,1]
df.rpos.allgen[,13] <- dir[,2]
dir <- str_split(df.rpos.allgen$data.dir, "-FILE", n=3, simplify = TRUE)
df.rpos.allgen[,11] <- dir[,1]

df.rpos.allgen[,12] <- as.numeric(df.rpos.allgen[,12])
df.rpos.allgen[,13] <- as.numeric(df.rpos.allgen[,13])

df.rpos10000 <- subset(df.rpos.allgen, Gen==10000)
df.rpos5000 <- subset(df.rpos.allgen, Gen==5000)
df.rpos5250 <- subset(df.rpos.allgen, Gen==5250)


ggplot(df.rpos5250, aes(x=V12, y=V13, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako")+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

ggplot(df.rpos5000, aes(x=V12, y=V13, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako")+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

ggplot(df.rpos10000, aes(x=V12, y=V13, z=ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = st_dev_abs)+
  # stat_summary_2d(breaks = c(-1.25,-0.75, -0.25, 0.25, 0.75, 1.25), fun = Angle_Mean)+
  scale_fill_viridis_c(option = "mako")+
  geom_point(aes(y=0, x=0))+
  labs(title = ("A/ Angle SD for different optimum \nphenotypes with no correlation of selection"), x="Expression gene A", y="Expression gene B")+
  labs(fill = "Ma sd")

ggplot(df.rpos5000, aes(x=V12, y=V13, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))

ggplot(df.rpos10000, aes(x=V12, y=V13, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))


ggplot(df.rpos5000, aes(x=V12, y=V13, z=Fitness))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("Fitness"))

ggplot(df.rpos10000, aes(x=V12, y=V13, z=Fitness))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("Fitness"))

ggplot(df.rpos5250, aes(x=V12, y=V13, z=Fitness))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("Fitness"))

ggplot(df.rpos5250, aes(x=V12, y=V13, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = TRUE)+
  scale_fill_viridis_c(option = "plasma")+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Free AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))
