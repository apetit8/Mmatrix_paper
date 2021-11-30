source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/ggplot2_custom_functions.R")
########################################################################################################################
sims.dirs.2 <- list.dirs(
  "../../simul/fig_4/init_netw", recursive = FALSE
)
#####################

df.all2 <- data.frame(NULL)
for (j in sims.dirs.2) {
  for(i in list.dirs(j, recursive = FALSE)){
    df <- df.raster.map(i, all.gen = TRUE)
    po <- str_split(i, "../../simul/fig_4/init_netw/", n=2, simplify = TRUE)
    pop <- str_split(po, "/angle_", n=2, simplify = TRUE)
    df[,9] <- sprintf("%s", pop[2,2])
    df[,10] <- sprintf("%s", pop[2,1])
    df.all2 <- rbind(df.all2, df)
  }
}

g3 <- ggplot(df.all2, aes(x=Gen, y=ang_M))+
  stat_summary(aes(color=V10, shape = V9), geom= "point", position = "identity", fun.data = Angle_Mean, size=3)+
  stat_summary(aes(x=10000, y=ang_S, shape = V9), geom= "point", position = "identity", fun.data = Angle_Mean, size=3)+
  labs(title = ("M correlation evolution for different initial correlation through generations"), x="Generations", y="M orientation")+
  labs(color = "Initial correlation", shape ="S orientation")+
  scale_shape_manual(labels = c("-pi/8", "-pi/2","pi/2", "0","3pi/8"), values = c(16,15,17,5,6))+
  scale_color_manual(labels = c("null", "negative","positive"), values = c("cornflowerblue", "orange","orchid4"))
g3 <- fracAx(p=g3, y=TRUE, x=FALSE, symbol = "pi")

pdf("../../figures/fig2-5_gen.pdf", width=9, height=5)
g3
dev.off()

########################################################################################################################
sims.dirs.3 <- list.dirs(
  "../../simul/fig_3/init_netw_full", recursive = FALSE
)
#####################

df.all3 <- data.frame(NULL)
for (j in sims.dirs.3) {
  for(i in list.dirs(j, recursive = FALSE)){
    df <- df.raster.map(i, all.gen = TRUE)
    po <- str_split(i, "../../simul/fig_3/init_netw_full/", n=2, simplify = TRUE)
    pop <- str_split(po, "/angle_", n=2, simplify = TRUE)
    df[,9] <- sprintf("%s", pop[2,2])
    df[,10] <- sprintf("%s", pop[2,1])
    df.all3 <- rbind(df.all3, df)
  }
}

for (i in unique(df.all3$V10)){
  print(i)
  ww <- subset(df.all3, V10 == i)
  delta_a <- 1 - mean((modulo.all(ww$ang_M-ww$ang_S, modulo=pi)^2)/((pi^2)/12))
  print(delta_a)
}