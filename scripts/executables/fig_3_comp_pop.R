source("../functions_R/network.R")
source("../functions_R/tools.R")
library(ggplot2)
library(ggExtra)
library(grid)
library(gridExtra)
library(ggstance)
library(gtable)
library(scales)
#This script trace 
#####################
sims.dirs <- c(
  "../../simul/fig_3/se_0","../../simul/fig_3/se_1",  "../../simul/fig_3/se_-1", "../../simul/fig_3/se1_ab_neg",
  "../../simul/fig_3/se0_16g",  "../../simul/fig_3/se0_p500", "../../simul/fig_3/se0_p10000"
)
of        <- "fig3_all"
modulo <- pi
#####################

df.all <- data.frame(NULL)
for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
  df <- df.topo.raw(sims.dir, network=FALSE)
  df[,9] <- modulo.all(df$ang_M-df$ang_S, modulo = modulo)
  mse <- round(mse.ms(df[,9]), 4)
  pop <- str_split(i, "../../simul/fig_3/", n=2, simplify = TRUE)
  df[,10] <- sprintf("%s", pop[,2])
  df.all <- rbind(df.all, df)
}
names(df.all)[names(df.all) == "V9"] <- "ang_diff"
names(df.all)[names(df.all) == "V10"] <- "MSE"

p <- ggplot()+
  labs(title = ("Absolute difference between M and S angle density"), y="Density", x="M S Angle difference")+
  # geom_histogram(data=df.all, aes(x=ang_diff, y=..density.., color=MSE), fill="white", binwidth = 0.05, alpha=0, position="identity")+
  geom_density(data=df.all, aes(abs(ang_diff), color=MSE, fill=MSE), alpha=0, lwd = 1, linetype = 1)

p0 <- ggplot()+
  labs(title = ("Absolute difference between M and S angle distribution"), y="Population", x="M S Angle difference")+
  geom_boxplot(data=df.all, aes(x=abs(ang_diff), y=MSE , fill=MSE), inherit.aes = FALSE)+
  scale_y_discrete(labels = NULL, breaks = NULL)

p1 <-  ggplot(data=df.all, aes(ang_S, ang_diff,col=MSE))+
  # coord_fixed(xlim = c(-pi/2,pi/2))+
  labs(title = ("Difference between M and S angle for different S angle"), x="S Angle", y="M S Angle difference")+
  geom_smooth(aes(fill=MSE), alpha = 0.1, method = "loess", span=0.15, level = 0.95)
p1 <- fracAx(p1, "pi")

p2 <- ggplot(data=df.all, aes(ang_S, ang_M,col=MSE))+
  coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
  geom_abline()+
  labs(title = ("M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  geom_smooth(aes(fill=MSE), alpha = 0.1, method = "loess", span=0.15, level = 0.95)+
  geom_point(aes(y=mean_ang_M, fill=MSE))
#       geom_line(aes(y=mean_ang_M))
p2 <- fracAx(p=p2, y=TRUE, "pi")

sims.dirs.wit <- list.dirs("../../simul/fig_2/round_s")
sims.dir.wit <- sims.dirs.wit[2:32]
diff.wit <- diff.ms.df(sims.dir.wit, modulo = modulo)

#Tab of MSE and statistic test between pop and witnesses
df.mse <- data.frame(ncol(4))
j <- 1
for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
  diff <- diff.ms.df(sims.dir, modulo = modulo)
  if (i==sims.dirs[1]){diff_se0 <- diff}
  mse <- round(mse.ms(diff),4)
  wx_se0 <- wilcox.test(((diff)^2),((diff_se0)^2))
  wx_wit <- wilcox.test(((diff)^2),((diff.wit)^2))
  pop <- str_split(i, "../../simul/fig_3/", n=2, simplify = TRUE)
  
  df.mse[j,1] <- pop[2]
  df.mse[j,2] <- mse
  df.mse[j,3] <-round(wx_se0$p.value,4)
  df.mse[j,4] <-round(wx_wit$p.value,4)
  j <- j +1
}

names(df.mse) <- c("Population","MSE","Wx test w/ se_0","Wx test w/ round S")

g <- tableGrob(as.matrix(df.mse),rows = NULL)

pdfname   <- print(sprintf("../../figures/%s_pop_comp.pdf", of))
# pdf(pdfname, width=16, height=10)
pdf(pdfname, width=16, height=12)
lay <- rbind(c(1,4),
             c(2,5),
             c(3,5))

grid.arrange(
  g,p,p0,p1,p2,
  ncol = 2,
  nrow = 3,
  widths = c(1,1),
  layout_matrix = lay,
  clip = FALSE
)

dev.off()
