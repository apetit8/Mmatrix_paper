source("../functions_R/network.R")
source("../functions_R/tools.R")
library(ggplot2)
library(ggExtra)
library(grid)
library(gridExtra)
library(ggstance)
#This script trace 
#####################
sims.dirs <- c(
         "../../simul/fig_3/se_0","../../simul/fig_3/se_1","../../simul/fig_3/se_-1", "../../simul/fig_3/se1_ab_neg", "../../simul/fig_2/round_s"
        # "../../simul/fig_2/c0_0_free","../../simul/fig_2/c0_0abba","../../simul/fig_2/c0_0adda", "../../simul/fig_2/c0_0ab","../../simul/fig_2/round_s"
         # "../../simul/fig_2/round_s", "../../simul/fig_2/round_s3", "../../simul/fig_2/witness_0w"
        )
of        <- "fig3_modolopi"
modulo <- pi
#####################


df.all <- data.frame(NULL)
for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
  df <- df.topo.raw(sims.dir, network=FALSE)
  df[,9] <- modulo.all(df$ang_M-df$ang_S, modulo = modulo)
  mse <- round(mse.ms(df[,9]), 4)
  pop <- str_split(i, "../../simul/fig_", n=2, simplify = TRUE)
  df[,10] <- sprintf("%s, MSE= %s", pop[,2], mse)
  df.all <- rbind(df.all, df)
}
names(df.all)[names(df.all) == "V9"] <- "ang_diff"
names(df.all)[names(df.all) == "V10"] <- "MSE"

pdfname   <- print(sprintf("../../figures/%s_pop_comp.pdf", of))
pdf(pdfname, width=8, height=17)

p <- ggplot()+
  labs(title = sprintf("%s, Absolute difference between M and S angle density", of), y="Density", x="M S Angle difference")+
  # geom_histogram(data=df.all, aes(x=ang_diff, y=..density.., color=MSE), fill="white", binwidth = 0.05, alpha=0, position="identity")+
  geom_density(data=df.all, aes(abs(ang_diff), color=MSE, fill=MSE), alpha=0.05, lwd = 1, linetype = 1)
  
p0 <- ggplot()+
  labs(title = sprintf("%s, Absolute difference between M and S angle distribution", of), y="Density", x="M S Angle difference")+
  geom_boxplot(data=df.all, aes(x=abs(ang_diff), y=MSE , fill=MSE), inherit.aes = FALSE)+
  scale_y_discrete(labels = NULL, breaks = NULL)

p1 <- ggplot(data=df.all, aes(ang_S, ang_diff,col=MSE))+
      # coord_fixed(xlim = c(-pi/2,pi/2))+
      labs(title = sprintf("%s, Difference between M and S angle for different S angle", of), x="S Angle", y="M S Angle difference")+
      geom_smooth(aes(fill=MSE), alpha = 0.1, method = "loess", span=0.15, level = 0.95)

p2 <- ggplot(data=df.all, aes(ang_S, ang_M,col=MSE))+
      coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
      geom_abline()+
      labs(title = sprintf("%s, M angle distribution and mean for different S angle", of), y="M Angle", x="S Angle")+
      geom_smooth(aes(fill=MSE), alpha = 0.1, method = "loess", span=0.15, level = 0.95)+
      geom_point(aes(y=mean_ang_M, fill=MSE))
#       geom_line(aes(y=mean_ang_M))

grid.newpage()
pushViewport(viewport(layout = grid.layout(6, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p, vp = vplayout(1, 1))
print(p0, vp = vplayout(2, 1))
print(p1, vp = vplayout(3, 1))  # key is to define vplayout
print(p2, vp = vplayout(4:5, 1))

dev.off()


# ggplot()+
#   labs(title = sprintf("%s, Difference between M and S angle", of), y="Density", x="M S Angle difference")+
#   # geom_histogram(data=df.all, aes(x=ang_diff, y=..density.., color=MSE), fill="white", binwidth = 0.05, alpha=0, position="identity")+
#   geom_density(data=df.all, aes(ang_diff, color=MSE, fill=MSE), alpha=0.05, lwd = 1, linetype = 1)+
#   geom_boxplot(data=df.all, aes(x=ang_diff, y=MSE , fill=MSE), inherit.aes = FALSE)+
#   stat_boxplot(geom = "vline", aes(xintercept = ..xlower..)) +
#   stat_boxplot(geom = "vline", aes(xintercept = ..xmiddle..)) +
#   stat_boxplot(geom = "vline", aes(xintercept = ..xupper..)) +
#   facet_grid(MSE ~ .) +
#   scale_fill_discrete()