source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/ggplot2_custom_functions.R")
library(png)
#This script trace 
#####################
sims.dirs <- c(
  "../../simul/fig_2/se_0","../../simul/fig_2/se_0_abba","../../simul/fig_2/se_0_o_abbacddc",
  "../../simul/fig_2/se_0_o_dacb" )
of        <- "fig2"
modulo <- pi
#####################

#Data
df.fig2 <- df.data(sims.dirs, pattern = "../../simul/fig_2/", variable="netw")

#Plot
p2 <- ggplot(data=df.fig2, aes(ang_S, ang_M,col=netw))+
  coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
  geom_abline()+
  labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  # geom_smooth(aes(fill=netw), alpha = 0, method = "AngleSmooth", span=0.9, level = 0.95)+
  geom_point(aes(y=mean_ang_M, fill=netw), show.legend = FALSE)+
  stat_summary(geom = "linerange", position = "identity", fun.data = st_dev_angle, size=3, alpha=0.2, show.legend = FALSE)
p2 <- fracAx(p=p2, y=TRUE, "pi")
p2 <- p2 + facet_wrap(vars(netw), ncol=2) + theme(strip.background = element_blank(), strip.text = element_blank())


#Tab
g_mean <- table.index(sims.dirs, pattern = "../../simul/fig_2/", ref = df.fig2, bymean = TRUE)
plot(g_mean)
g_sd <- table.index(sims.dirs, pattern = "../../simul/fig_2/", ref = df.fig2, bymean = FALSE)
plot(g_sd)


#PNG
png.netw1 = readPNG('../../simul/fig_2/network1.png')
png.netw2 = readPNG('../../simul/fig_2/network2.png')
png.netw3 = readPNG('../../simul/fig_2/network3.png')
png.netw4 = readPNG('../../simul/fig_2/network4.png')

###############PDF
pdfname   <- print(sprintf("../../figures/%s.pdf", of))
pdf(pdfname, width=10, height=7)

grid.arrange(
  p2,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)

grid.raster(png.netw1, x=.1, y=0.75, width=0.14)
grid.raster(png.netw2, x=.93, y=0.75, width=0.14)
grid.raster(png.netw3, x=.1, y=0.3, width=0.14)
grid.raster(png.netw4, x=.93, y=0.3, width=0.14)

dev.off()
grid.arrange(g)
