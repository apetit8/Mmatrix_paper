source("../functions_R/All_functions.R")
########################################################################################################################
sims.dir.mutrate <- list.dirs(
  "../../simul/supp_data/mutrate", recursive = FALSE
)
of        <- "supp_data3"
modulo <- pi
#####################

df.mutrate <- df.data(sims.dir.mutrate, pattern = "../../simul/supp_data/mutrate/se0_mutrate_", variable="mutrate")

pm <- ggplot(data=df.mutrate, aes(ang_S, ang_M,col=mutrate))+
  coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
  geom_abline()+
  labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  # geom_smooth(aes(fill=mutrate), alpha = 0, method = "AngleSmooth", span=0.9, level = 0.95)+
  geom_point(aes(y=mean_ang_M, fill=mutrate), show.legend = FALSE)+
  stat_summary(geom = "linerange", position = "identity", fun.data = st_dev_angle, size=3, alpha=0.2, show.legend = FALSE)
pm <- fracAx(p=pm, y=TRUE, "pi")
pm <- pm + facet_wrap(vars(mutrate), ncol=2)
pm

##Tab of mean, standard deviation and alignment index
g_meanmr <- table.index(sims.dir.mutrate, pattern = "../../simul/supp_data/mutrate/se0_mutrate_", ref = df.mutrate, bymean = TRUE)
g_sdmr <- table.index(sims.dir.mutrate, pattern = "../../simul/supp_data/mutrate/se0_mutrate_", ref = df.mutrate, bymean = FALSE)

source("../functions_R/test.R")
pmutrate <- ggplot(data=df.mutrate, aes(ang_S, ang_M,col=mutrate))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline()+
  geom_abline(intercept=pi)+
  geom_abline(intercept=-pi)+
  geom_point(aes(y=ang_M, fill=mutrate), alpha=0.1, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, fill=mutrate), alpha=0.1, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, fill=mutrate), alpha=0.1, show.legend = FALSE)+
  labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  # geom_smooth(aes(fill=mutrate), alpha = 0, method = "AngleSmooth", span=0.9, level = 0.95)+
  geom_point(aes(y=mean_ang_M, fill=mutrate), show.legend = FALSE)
# stat_summary(geom = "linerange", position = "identity", fun.data = st_dev_angle, size=3, alpha=0.2, show.legend = FALSE)
pmutrate <- fracAx(p=pmutrate, y=TRUE, "pi")
pmutrate <- pmutrate + facet_wrap(vars(mutrate), ncol=2) + theme(strip.background = element_blank())
pmutrate

###############PDF
pdfname   <- print(sprintf("../../figures/%s_mutrate.pdf", of))
cairo_pdf(pdfname, width=10, height=7)

lay <- rbind(c(1,3),
             c(1,2))

grid.arrange(
  pmutrate, g_meanmr, g_sdmr,
  ncol = 2,
  nrow = 2,
  widths = c(1,0.75),
  layout_matrix = lay,
  clip = FALSE
)

dev.off()
grid.arrange(g)



# df.mutgen <- data.frame(NULL)
# for (j in sims.dir.mutrate) {
#   for(i in list.dirs(j, recursive = FALSE)){
#     df <- df.raster.map(i, all.gen = TRUE)
#     po <- str_split(i, "../../simul/supp_data/mutrate/se0_mutrate_", n=2, simplify = TRUE)
#     pop <- str_split(po, "/angle_", n=2, simplify = TRUE)
#     df[,9] <- sprintf("%s", pop[2,2])
#     df[,10] <- sprintf("%s", pop[2,1])
#     df.mutgen <- rbind(df.mutgen, df)
#   }
# }
#
# dfmg <- subset(df.mutgen, V9 == "-0.7.par")
# g3 <- ggplot(dfmg, aes(x=Gen, y=ang_M))+
#   stat_summary(aes(color=V10), geom= "point", position = "identity", fun.data = Angle_Mean, size=3)
#   # stat_summary(aes(x=10000, y=ang_S, shape = V9), geom= "point", position = "identity", fun.data = Angle_Mean, size=3)+
#   # labs(title = ("M correlation evolution for different initial correlation through generations"), x="Generations", y="M orientation")+
#   # labs(color = "Initial correlation", shape ="S orientation")
#   # scale_shape_manual(labels = c("-pi/8", "-pi/2","pi/2", "0","3pi/8"), values = c(16,15,17,5,6))+
#   # scale_color_manual(labels = c("null", "negative","positive"), values = c("cornflowerblue", "orange","orchid4"))
# g3 <- fracAx(p=g3, y=TRUE, x=FALSE, symbol = "pi")
# g3