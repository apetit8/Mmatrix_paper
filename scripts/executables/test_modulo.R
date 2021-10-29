modulo.all <- function(angle, modulo=pi) {
  angle <- angle %% modulo
  ifelse(angle > modulo/2, angle - modulo, angle)
}

mean.angle.2pi <- function(angles) {
  atan2(mean(sin(angles)), mean(cos(angles)))
}

mean.angle.pi <- function(angles) {
  modulo.all(mean.angle.2pi(2*(angles %% pi))/2, pi)
}



df.all <- data.frame(NULL)
for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:30]
  df <- df.topo.raw(sims.dir, network=FALSE)
  df[,8] <- modulo.all(df$ang_M-df$ang_S, modulo=pi/2)
  mse <- round(mse.ms(df[,8]), 4)
  pop <- str_split(i, "../../simul/fig_", n=2, simplify = TRUE)
  df[,9] <- sprintf("%s, MSE= %s", pop[,2], mse)
  df.all <- rbind(df.all, df)
}
names(df.all)[names(df.all) == "V8"] <- "ang_diff"
names(df.all)[names(df.all) == "V9"] <- "MSE"


mean.angle.2pi()

ggplot()+
  labs(title = sprintf("%s, Difference between M and S angle", of), y="Density", x="M S Angle difference")+
  # geom_histogram(data=df.all, aes(x=ang_diff, y=..density.., color=MSE), fill="white", binwidth = 0.05, alpha=0, position="identity")+
  geom_density(data=df.all, aes(ang_diff, color=MSE, fill=MSE), alpha=0.05, lwd = 1, linetype = 1)

ggplot()+
  labs(title = sprintf("%s, Difference between M and S angle", of), y="Density", x="M S Angle difference")+
  geom_boxplot(data=df.all, aes(x=ang_diff, y=MSE , fill=MSE), inherit.aes = FALSE)+
  geom_density(data=df.all, aes(ang_diff_arnd, color=MSE, fill=MSE), alpha=0.05, lwd = 1, linetype = 3)+
  scale_y_discrete(labels = NULL, breaks = NULL)


ggplot(data=df.all, aes(ang_S, ang_diff,col=MSE))+
  # coord_fixed(xlim = c(-pi/2,pi/2))+
  labs(title = sprintf("%s, Difference between M and S angle for different S angle", of), x="S Angle", y="M S Angle difference")+
  geom_smooth(aes(fill=MSE), alpha = 0.3, span=0.1, level = 0.95)+
  geom_smooth(aes(y=abs(ang_diff_arnd), fill=MSE), alpha = 0.3, span=0.3, level = 0.95, linetype = 3)
  
