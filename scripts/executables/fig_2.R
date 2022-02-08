source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs <- c(
  "../../simul/fig_2/se_0","../../simul/fig_2/se_0_abba","../../simul/fig_2/se_0_o_abbacddc",
  "../../simul/fig_2/se_0_o_dacb" )
modulo <- pi
#####################

#Data
df.fig2 <- df.data(sims.dirs, pattern = "../../simul/fig_2/", variable="netw")

pfig2 <- ggplot(data=df.fig2, aes(ang_S, ang_M,col=netw))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  #previously : alpha of 0.1
  geom_point(aes(y=ang_M, fill=1-ecc_M), alpha=0.12, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, fill=1-ecc_M), alpha=0.12, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, fill=1-ecc_M), alpha=0.12, show.legend = FALSE)+
  # geom_point(aes(y=mean_ang_M, fill=netw), show.legend = FALSE)+
  labs(y="M direction", x="S direction")
pfig2 <- pfig2 + scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                                  labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                       labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig2 <- pfig2 + facet_wrap(vars(netw), ncol=2) + theme(strip.background = element_blank(), strip.text = element_blank())
pfig2

#Tab
# g_meanfig2 <- table.index(sims.dirs, pattern = "../../simul/fig_2/", ref = df.fig2, bymean = TRUE)
# plot(g_meanfig2)
g_sdfig2 <- table.index(sims.dirs, pattern = "../../simul/fig_2/", ref = df.fig2, bymean = FALSE)
plot(g_sdfig2)
# 
# df.4 <- subset(df.fig2, netw == "se_0_o_dacb")
# 1 - (1/(((pi^2)/12)*nrow(df.4)))*sum(df.4$ang_diff^2)
# mean(df.4$ang_diff^2)
# (pi^2)/12
# #Mean square distance
# sum(df.4$ang_diff^2)/nrow(df.4) #Why th is it not closer to (pi^2)/12 ??! Like the other one where points where more widely distributed ?

#PNG
png.netw1 = readPNG('../../simul/fig_2/network1.png')
png.netw2 = readPNG('../../simul/fig_2/network2.png')
png.netw3 = readPNG('../../simul/fig_2/network3.png')
png.netw4 = readPNG('../../simul/fig_2/network4.png')

###############PDF
cairo_pdf("../../figures/fig_2.pdf", width=10, height=8)

grid.arrange(
  pfig2,
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


cairo_pdf("../../figures/fig_2_small.pdf", width=7.5, height=6)

grid.arrange(
  pfig2,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
grid.raster(png.netw1, x=.1, y=0.75, width=0.14)
grid.raster(png.netw2, x=.92, y=0.75, width=0.14)
grid.raster(png.netw3, x=.1, y=0.3, width=0.14)
grid.raster(png.netw4, x=.92, y=0.3, width=0.14)

dev.off()




# 
# 
# pfig2 <- ggplot(data=df.fig2, aes(ang_S, ang_M,col=ecc_M))+
#   coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
#   geom_abline()+
#   geom_abline(intercept=pi)+
#   geom_abline(intercept=-pi)+
#   #previously : alpha of 0.1
#   geom_point(aes(y=ang_M, col=ecc_M), alpha=0.3)+
#   geom_point(aes(y=ang_M_mpi, col=ecc_M), alpha=0.3)+
#   geom_point(aes(y=ang_M_ppi, col=ecc_M), alpha=0.3)+
#   # geom_point(aes(y=mean_ang_M, fill=netw), show.legend = FALSE)+
#   labs(title = ("A/ M direction distribution and mean for different S direction"), y="M direction", x="S direction")
# pfig2 <- pfig2 + scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
#                                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
#   scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
#                      labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
# pfig2 <- pfig2 + facet_wrap(vars(netw), ncol=2) + theme(strip.background = element_blank(), strip.text = element_blank())
# pfig2
