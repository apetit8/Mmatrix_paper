source("../functions_R/All_functions.R")
########################################################################################################################
sims.dirs <- "../../simul/fig_3_bis/no_ab"
of        <- "fig3"
modulo <- pi
#####################

df.rnul <- df.opt.map2(sims.dirs, modulo=modulo)

g1 <- ggplot(df.rnul, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = FALSE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  geom_point(aes(y=0, x=0))+
  labs(title =  ("A/ Without initial AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))
g1

######################
sims.dirs.5 <- "../../simul/fig_3_bis/neg_ab"
modulo <- pi
#####################

df.rneg <- df.opt.map2(sims.dirs.5, gen=FALSE)

g3 <-   ggplot(df.rneg, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean")+
  geom_point(aes(col=ang_M))+
  # scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  scale_fill_viridis_c(option = "plasma")+
  labs(title =  ("C/ Negative initial AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))
g3



df.rneg22 <- df.opt.map2(sims.dirs.5, gen=TRUE)



# ggplot(df.rneg.4500, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
#   coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
#   stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean")+
#   # scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
#   labs(title =  ("C/ Negative initial AB regulation"), x="Expression gene A", y="Expression gene B")+
#   labs(fill = expression("\u03BE\u03B1"))

pdfname   <- print("../../figures/neg_corr_AB_fixed2.pdf")
pdf(pdfname, width=5, height=5)

for (i in c(0,500,1000,2000,4000,5000,5250,6000,8000,10000,13000,15000)) {
  
  df.rneg.4500 <- subset(df.rneg22, Gen == i)
  pp <-  ggplot(df.rneg.4500, aes(x=P_mean_A, y=P_mean_B, col=ang_M))+
    geom_point()+
    labs(title =  i)
  print(pp)
}

dev.off()

######################
sims.dirs.6 <- "../../simul/fig_3_bis/pos_ab"
modulo <- pi
#####################

df.rpos<- df.opt.map(sims.dirs.6)

g4 <- ggplot(df.rpos, aes(x=P_mean_A, y=P_mean_B, z=sq_dist))+
  coord_fixed(ratio = 1, xlim = c(-1.125,1.125), ylim = c(-1.125,1.125), expand = TRUE, clip = "on")+
  stat_summary_2d(breaks = seq(-1.125, 1.125, by = 0.25), fun = "mean", show.legend = FALSE)+
  scale_fill_viridis_c(option = "plasma", limits=c(-0.1, 1))+
  labs(title =  ("B/ Positive initial AB regulation"), x="Expression gene A", y="Expression gene B")+
  labs(fill = expression("\u03BE\u03B1"))
g4

#################################
# PDF
#################################


