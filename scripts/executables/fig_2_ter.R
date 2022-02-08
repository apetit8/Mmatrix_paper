source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs.fig2 <-list.dirs(
  "../../simul/fig_2_ter", recursive = FALSE
)

# sims.dirs.fig2 <-c("../../simul/fig_2_bis/d1","../../simul/fig_2_bis/d3", #"../../simul/fig_2_bis/d1",
# "../../simul/fig_2_bis/d4","../../simul/fig_2_bis/d5")
modulo <- pi
#####################

#Data
df.fig2c <- df.data(sims.dirs.fig2, pattern = "../../simul/fig_2_ter/g", variable="netw")
g_sdfig2b <- table.index(sims.dirs.fig2, pattern = "../../simul/fig_2_ter/g", ref = df.fig2c, bymean = FALSE, asgrob=FALSE)
g_sdfig2b[,5] <- c( mean(subset(df.fig2c, netw == "12")$ecc_M),
                    mean(subset(df.fig2c, netw == "16")$ecc_M),
                    mean(subset(df.fig2c, netw == "4")$ecc_M),
                    mean(subset(df.fig2c, netw == "6")$ecc_M),
                    mean(subset(df.fig2c, netw == "8")$ecc_M))
names(g_sdfig2b) <- c("Population","sq_rho","standard deviation","Xi_alpha", "Eccentricity")
g_sdfig2b[,1] <-c(12,16,4,6,8)

dir <- str_split(df.fig2c$data.dir, "simul", n=3, simplify = TRUE)
df.fig2c[,11] <- dir[,2]
df_angevolv <- data.frame()
j <-1
for (u in as.list(df.fig2c[,11])){
  dist <-  1- mse.ms(subset(df.fig2c, V11 == u)$ang_diff)/((pi^2)/12)
  df_angevolv[j,1] <- dist
  j <-j+1
}
df.fig2c[,12] <- df_angevolv

#Plot
pfig2b <- ggplot(data=g_sdfig2b, aes(Population, Xi_alpha))+
  geom_col()+
  labs(y="\u03BE\u03B1", x="Network Size")
pfig2b

ggplot(data=df.fig2c, aes(as.numeric(netw), sqrt(ang_diff^2) ) )+
  geom_point()+
  labs(y="\u03BE\u03B1", x="Network Size")


###############PDF
cairo_pdf("../../figures/fig_2_bis.pdf", width=6, height=3)
grid.arrange(
  pfig2b,
  ncol = 1,
  nrow = 1,
  clip = FALSE
)

grid.raster(png.dist.netw, x=0.75, y=0.5, width=0.5)

dev.off()


###########

dir <- str_split(df.fig2$data.dir, "simul", n=3, simplify = TRUE)
df.fig2[,11] <- dir[,2]
df_angevolv <- data.frame()
j <-1
for (u in as.list(df.fig2[,11])){
  dist <-  1- mse.ms(subset(df.fig2, V11 == u)$ang_diff)/((pi^2)/12)
  df_angevolv[j,1] <- dist
  j <-j+1
}
df.fig2[,12] <- df_angevolv


g_ntwsize <- ggplot(data=df.fig2c, aes(as.numeric(netw), V1) )+
  geom_point(aes(col=ang_S ), show.legend = FALSE)+
  stat_summary(fun=mean, colour="red", geom="line")+
  scale_color_gradient2(low="blue", mid="yellow",
                        high="blue", space ="Lab" )+
  labs(y="\u03BE\u03B1", x="Network Size")+
  theme_light()

g_conect <- ggplot(data=df.fig2, aes(netw, V1) )+
  geom_point(aes(col=ang_S ))+
  stat_summary(fun=mean, colour="red", geom="line")+
  scale_color_gradient2(low="blue", mid="yellow",
                        high="blue", space ="Lab" )+
  labs(y="\u03BE\u03B1", x="A B connectivity")+
  theme_light()+
  theme(plot.margin = unit(c(0.25,6, 0.2, 0.1), "cm"))
  
arrow = readPNG('../../templates/fig_2_bis/arrow.png')


cairo_pdf("../../figures/fig_2_ter2.pdf", width=9, height=3)
grid.arrange(
  g_ntwsize, g_conect,
  ncol = 2,
  nrow = 1,
  widths = c(0.5, 1),
  clip = FALSE
)
# grid.raster(png.dist.netw, x=0.84, y=0.5, width=0.33)
grid.raster(png.dist.netw, x=0.87, y=0.5, width=0.25)
grid.raster(arrow, x=0.255, y=0.07, width=0.28)
dev.off()





#################Circular plot ?
df.fig2[,13] <-  cos(df.fig2[,3])^2
df.fig2[,14] <- sign(df.fig2[,3])*sqrt(abs(df.fig2[,12]^2-df.fig2[,13]))


ggplot(data=df.fig2, aes(V1, V14) )+
  geom_point(aes(col=netw ))+
  theme_light()

  
  ggplot(data=df.fig2, aes(ang_S, V1) )+
    geom_point(aes(col=netw ))+
    geom_smooth(aes(col=netw ), alpha=0)+
    coord_polar(start = -pi/2, direction = 1)+
    scale_x_continuous(limits=c(-pi,pi),breaks=seq(-pi/2, pi/2, pi/4),labels=c("-\u03c0/2", "-\u03c0/4", "0", "\u03c0/4", "\u03c0/2"))+
    labs(x="\u03BE\u03B1", y="S direction")+
    theme_bw()
