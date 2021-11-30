source("../functions_R/figure_tools.R")

#####################
sims.dirs <- list.dirs("../../simul/fig_1")
#angle = c(-0.2, 1.2)
angle = c(1,2,3) #c(-1.4, -1.3, -1.2 , -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4)
of        <- "fig1"
modulo <- pi
#####################

df.m.s1 <- data.frame()
for (i in angle) {

  sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle__%s.par", i), full.names=TRUE)
  df <- df.raster.map(sims.dir, all.gen=TRUE)
  pop <- str_split(df$data.dir, "../../simul/fig_1/", n=2, simplify = TRUE)
  df[,9] <- sprintf("%s", pop[,2])
  model <- str_split(pop[,2], "/angle__", n=2, simplify = TRUE)
  df[,10] <- sprintf("%s", model[,1])
  df.m.s1 <- rbind(df.m.s1, df) 
  
}


df.m.s <- data.frame()
for (i in angle) {
  
  sims.dir  <- list.files(path=sims.dirs, pattern=sprintf("angle__%s.par", i), full.names=TRUE)
  df <- df.raster.map(sims.dir, all.gen=FALSE)
  pop <- str_split(df$data.dir, "../../simul/fig_1/", n=2, simplify = TRUE)
  df[,9] <- sprintf("%s", pop[,2])
  model <- str_split(pop[,2], "/angle__", n=2, simplify = TRUE)
  df[,10] <- sprintf("%s", model[,1])
  df.m.s <- rbind(df.m.s, df) 
  
}
#Calcul des evolvability for all angles at the same time

df.m <- subset(df.m.s, V10 == "m_2000")
df.w <- subset(df.m.s, V10 == "w2")

m_Eva <-round( 1- mean((modulo.all((df.m$mean_ang_m - df.m$ang_S) ,modulo=pi)^2)/((pi^2)/12)),3)
w_Eva <-round( 1- mean((modulo.all((df.w$mean_ang_m - df.m$ang_S) ,modulo=pi)^2)/((pi^2)/12)),3)

sprintf("\u03BE\u03B1 Wagner = %s", w_Eva)
sprintf("\u03BE\u03B1 Multilinear = %s", m_Eva)

gg <- ggplot(df.m.s1, aes(x=Gen, y=ang_M))+
  geom_smooth(aes(color=V10, linetype = as.character(round(ang_S, 4))), alpha = 0.1, method = AngleSmooth, span=0.15, level = 1)+
  geom_point(aes(x=10000, y=round(ang_S, 4)), size=2)+
  labs(title = ("M correlation evolution for different initial correlation through generations"), x="Generations", y="M orientation")+
  labs(linetype = "S direction", color ="Model")+
  scale_linetype_manual(labels = c("0","pi/8", "pi/2"), values = c(3,2,1))+
  scale_color_manual(labels = c("Multilinear", "Wagner"), values = c("cornflowerblue", "orange"))
gg <- fracAx(p=gg, y=TRUE, x=FALSE, symbol = "pi")





