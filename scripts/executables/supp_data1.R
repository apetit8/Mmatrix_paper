source("../functions_R/All_functions.R")
########################################################################################################################
sims.dirs <- list.dirs(
  "../../simul/supp_data/genes_nbr", recursive = FALSE
)
of        <- "supp_data"
modulo <- pi
#####################

df.genes_nbr <- df.data(sims.dirs, pattern = "../../simul/supp_data/genes_nbr/se0_", variable="genenbr")

pg <- ggplot(data=df.genes_nbr, aes(ang_S, ang_M,col=genenbr))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline()+
  geom_abline(intercept=pi)+
  geom_abline(intercept=-pi)+
  geom_point(aes(y=ang_M, fill=genenbr), alpha=0.1, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, fill=genenbr), alpha=0.1, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, fill=genenbr), alpha=0.1, show.legend = FALSE)+
  labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  geom_point(aes(y=mean_ang_M, fill=genenbr), show.legend = FALSE)
pg <- fracAx(p=pg, y=TRUE, "pi")
pg <- pg + facet_wrap(vars(genenbr), ncol=2) + theme(strip.background = element_blank())
pg

g_mean_gnb <- table.index(sims.dirs, pattern = "../../simul/supp_data/genes_nbr/se0_", ref = df.genes_nbr, bymean = TRUE)
plot(g_mean_gnb)
g_sd_gnb <- table.index(sims.dirs, pattern = "../../simul/supp_data/genes_nbr/se0_", ref = df.genes_nbr, bymean = FALSE)
plot(g_sd_gnb)

###############PDF
pdfname   <- print(sprintf("../../figures/%s_nbrgenes.pdf", of))
cairo_pdf(pdfname, width=10, height=7)

lay <- rbind(c(1,3),
             c(1,2))

grid.arrange(
  pg, g_mean_gnb, g_sd_gnb,
  ncol = 2,
  nrow = 2,
  widths = c(1,0.75),
  layout_matrix = lay,
  clip = FALSE
)

dev.off()

########################################################################################################################
sims.dir.pop <- list.dirs(
  "../../simul/supp_data/popsize", recursive = FALSE
)
of        <- "supp_data"
#####################

df.popsize <- df.data(sims.dir.pop, pattern = "../../simul/supp_data/popsize/se0_p", variable="pops")

pp <- ggplot(data=df.popsize, aes(ang_S, ang_M,col=pops))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline()+
  geom_abline(intercept=pi)+
  geom_abline(intercept=-pi)+
  geom_point(aes(y=ang_M, fill=pops), alpha=0.1, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, fill=pops), alpha=0.1, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, fill=pops), alpha=0.1, show.legend = FALSE)+
  labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  geom_point(aes(y=mean_ang_M, fill=pops), show.legend = FALSE)
pp <- fracAx(p=pp, y=TRUE, "pi")
pp <- pp + facet_wrap(vars(pops), ncol=2) + theme(strip.background = element_blank())
pp

##Tab of mean, standard deviation and alignment index
g_mean_ps <- table.index(sims.dir.pop, pattern = "../../simul/supp_data/popsize/se0_p", ref = df.popsize, bymean = TRUE)
g_sd_ps <- table.index(sims.dir.pop, pattern = "../../simul/supp_data/popsize/se0_p", ref = df.popsize, bymean = FALSE)

###############PDF
pdfname   <- print(sprintf("../../figures/%s_popsize.pdf", of))
cairo_pdf(pdfname, width=10, height=7)

lay <- rbind(c(1,3),
             c(1,2))

grid.arrange(
  pp, g_mean_ps, g_sd_ps,
  ncol = 2,
  nrow = 2,
  widths = c(1,0.75),
  layout_matrix = lay,
  clip = FALSE
)

dev.off()

