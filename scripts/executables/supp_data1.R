source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/ggplot2_custom_functions.R")
########################################################################################################################
sims.dirs <- list.dirs(
  "../../simul/supp_data/genes_nbr", recursive = FALSE
)
of        <- "supp_data1"
modulo <- pi
#####################

df.genes_nbr <- df.data(sims.dirs, pattern = "../../simul/supp_data/genes_nbr/se0_", variable="genenbr")

pg <- ggplot(data=df.genes_nbr, aes(ang_S, ang_M,col=genenbr))+
  coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
  geom_abline()+
  labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  # geom_smooth(aes(fill=genenbr), alpha = 0, method = "AngleSmooth", span=0.9, level = 0.95)+
  geom_point(aes(y=mean_ang_M, fill=genenbr), show.legend = FALSE)+
  stat_summary(geom = "linerange", position = "identity", fun.data = st_dev_angle, size=3, alpha=0.2, show.legend = FALSE)
pg <- fracAx(p=pg, y=TRUE, "pi")
pg <- pg + facet_wrap(vars(genenbr), ncol=2)


g_mean <- table.index(sims.dirs, pattern = "../../simul/supp_data/genes_nbr/se0_", ref = df.genes_nbr, bymean = TRUE)
plot(g_mean)
g_sd <- table.index(sims.dirs, pattern = "../../simul/supp_data/genes_nbr/se0_", ref = df.genes_nbr, bymean = FALSE)
plot(g_sd)

########################################################################################################################
sims.dir.pop <- list.dirs(
  "../../simul/supp_data/popsize", recursive = FALSE
)
of        <- "supp_data2"
#####################

df.popsize <- df.data(sims.dir.pop, pattern = "../../simul/supp_data/popsize/se0_p", variable="pops")

pp <- ggplot(data=df.popsize, aes(ang_S, ang_M,col=pops))+
  coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
  geom_abline()+
  labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  # geom_smooth(aes(fill=pops), alpha = 0, method = "AngleSmooth", span=0.9, level = 0.95)+
  geom_point(aes(y=mean_ang_M, fill=pops), show.legend = FALSE)+
  stat_summary(geom = "linerange", position = "identity", fun.data = st_dev_angle, size=3, alpha=0.2, show.legend = FALSE)
pp <- fracAx(p=pp, y=TRUE, "pi")
pp <- pp + facet_wrap(vars(pops), ncol=2)

##Tab of mean, standard deviation and alignment index
g_mean1 <- table.index(sims.dir.pop, pattern = "../../simul/supp_data/popsize/se0_p", ref = df.popsize, bymean = TRUE)
g_sd1 <- table.index(sims.dir.pop, pattern = "../../simul/supp_data/popsize/se0_p", ref = df.popsize, bymean = FALSE)

###############PDF
pdfname   <- print(sprintf("../../figures/%s_popsize.pdf", of))
cairo_pdf(pdfname, width=10, height=7)

lay <- rbind(c(1,3),
             c(1,2))

grid.arrange(
  pp, g_mean1, g_sd1,
  ncol = 2,
  nrow = 2,
  widths = c(1,0.75),
  layout_matrix = lay,
  clip = FALSE
)

dev.off()
grid.arrange(g)

########################################################################################################################
sims.dir.mutrate <- list.dirs(
  "../../simul/supp_data/mutrate", recursive = FALSE
)
of        <- "supp_data3"
modulo <- pi
#####################

df.mutrate <- df.data(sims.dir.mutrate, pattern = "../../simul/supp_data/mutrate/se0_mutrate_", variable="mutrate")

pp <- ggplot(data=df.mutrate, aes(ang_S, ang_M,col=mutrate))+
  coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
  geom_abline()+
  labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  # geom_smooth(aes(fill=mutrate), alpha = 0, method = "AngleSmooth", span=0.9, level = 0.95)+
  geom_point(aes(y=mean_ang_M, fill=mutrate), show.legend = FALSE)+
  stat_summary(geom = "linerange", position = "identity", fun.data = st_dev_angle, size=3, alpha=0.2, show.legend = FALSE)
pp <- fracAx(p=pp, y=TRUE, "pi")
pp <- pp + facet_wrap(vars(mutrate), ncol=2)

##Tab of mean, standard deviation and alignment index
g_mean2 <- table.index(sims.dir.mutrate, pattern = "../../simul/supp_data/mutrate/se0_mutrate_", ref = df.mutrate, bymean = TRUE)
g_sd2 <- table.index(sims.dir.mutrate, pattern = "../../simul/supp_data/mutrate/se0_mutrate_", ref = df.mutrate, bymean = FALSE)
