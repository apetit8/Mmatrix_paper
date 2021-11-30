source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/ggplot2_custom_functions.R")
library(ggplot2)
library(ggExtra)
library(grid)
library(gridExtra)
library(ggstance)
library(gtable)
library(scales)
library(png)
#This script trace 
#####################
sims.dirs <- c(
         "../../simul/fig_2/se_0","../../simul/fig_2/se_0_abba","../../simul/fig_2/se_0_o_abbacddc",
         "../../simul/fig_2/se_0_o_dacb"
         # "../../simul/fig_2/se_0_ab", 
         # "../../simul/fig_2/se_0_ba",  "../../simul/fig_2/se_0_adda", "../../simul/fig_2/se_0_ad",
         # "../../simul/fig_2/se_0_da",
         # "../../simul/fig_2/round_s", "../../simul/fig_2/round_s3", "../../simul/fig_2/witness_0w"
        )
of        <- "fig2"
modulo <- pi
#####################

df.all <- data.frame(NULL)
for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
  df <- df.topo.raw(sims.dir, network=FALSE)
  df[,9] <- modulo.all(df$ang_M-df$ang_S, modulo = modulo)
  sq_rho <- round(mse.ms(df[,9]), 4)
  pop <- str_split(i, "../../simul/fig_2/", n=2, simplify = TRUE)
  df[,10] <- sprintf("%s", pop[,2])
  df.all <- rbind(df.all, df)
}
names(df.all)[names(df.all) == "V9"] <- "ang_diff"
names(df.all)[names(df.all) == "V10"] <- "netw"

p <- ggplot()+
  labs(title = ("B/ Difference between M and S angle density"), y="Density", x="M S Angle difference")+
  # geom_histogram(data=df.all, aes(x=ang_diff, y=..density.., color=netw), fill="white", binwidth = 0.05, alpha=0, position="identity")+
  geom_density(data=df.all, aes(ang_diff, color=netw, fill=netw), alpha=0, lwd = 1, linetype = 1, show.legend = FALSE)
p <- fracAx(p,y=FALSE,"pi")
  
# p0 <- ggplot()+
#   labs(title = ("Difference between M and S angle distribution"), y="Population", x="M S Angle difference")+
#   geom_boxplot(data=df.all, aes(x=ang_diff, y=netw , fill=netw), inherit.aes = FALSE)+
#   scale_y_discrete(labels = NULL, breaks = NULL)
# 
# p1 <-  ggplot(data=df.all, aes(ang_S, ang_diff,col=netw))+
#       # coord_fixed(xlim = c(-pi/2,pi/2))+
#       labs(title = ("Difference between M and S angle for different S angle"), x="S Angle", y="M S Angle difference")+
#       geom_smooth(aes(fill=netw), alpha = 0.1, method = "loess", span=0.15, level = 0.95)
# p1 <- fracAx(p1, "pi")

p2 <- ggplot(data=df.all, aes(ang_S, ang_M,col=netw))+
      coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
      geom_abline()+
      labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
      # geom_smooth(aes(fill=netw), alpha = 0, method = "AngleSmooth", span=0.9, level = 0.95)+
      geom_point(aes(y=mean_ang_M, fill=netw), show.legend = FALSE)+
      stat_summary(geom = "linerange", position = "identity", fun.data = st_dev_angle, size=3, alpha=0.2, show.legend = FALSE)
p2 <- fracAx(p=p2, y=TRUE, "pi")
p2 <- p2 + facet_wrap(vars(netw), ncol=2) + theme(strip.background = element_blank(), strip.text = element_blank())


##########Tab
df.sq_rho <- data.frame(ncol(4))
j <- 1
modulo <- pi
for (i in sims.dirs) {
  sims.dir <- list.dirs(i, recursive = FALSE)
  diff <- diff.ms.df(sims.dir, modulo = modulo)
  sq_rho <-round( mse.ms(diff), 4)
  delta_diff <-mean(mse.ms(diff)/(((pi^2)/12))) #distance if M angle is constant

  #mean square distance
  sdall <- data.frame()
  sddf <- filter(df.all, grepl(i,sprintf("%s/", data.dir)))
  yi <- unique(round(sddf$ang_S, 3))
  for (l in yi) {
    ww <- subset(sddf, round(ang_S, 3) == l)
    sdy <- as.data.frame(sqrt( modulo.all((ww$ang_M - mean.angle.pi(ww$ang_M)))^2))
    sdall <- rbind(sdall, sdy)
  }
  sd <- mean(sdall[,1]) #square deviation
  
  pop <- str_split(i, "../../simul/fig_2/", n=2, simplify = TRUE)
  df.sq_rho[j,1] <- pop[2]
  df.sq_rho[j,2] <- sq_rho
  df.sq_rho[j,3] <- round(sd, 4)
  df.sq_rho[j,4] <- round(1 - delta_diff, 4)
  j <- j +1
}
names(df.sq_rho) <- c("Population","sq_rho","standard deviation","Delta_a")
g <- tableGrob(as.matrix(df.sq_rho),rows = NULL)

png.netw = readPNG('../../simul/fig_2/Figure2_network.png')


pdfname   <- print(sprintf("../../figures/%s_pop_comp.pdf", of))
pdf(pdfname, width=10, height=6)

lay <- rbind(c(3,NA),
             c(3,2))
grid.arrange(
  p,p2,
  ncol = 2,
  nrow = 2,
  widths = c(1.2,1),
  layout_matrix = lay,
  clip = FALSE
)

grid.raster(png.netw, x=.8, y=0.75, width=0.3)

dev.off()
grid.arrange(g)
