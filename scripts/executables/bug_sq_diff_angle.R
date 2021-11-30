#####################
sims.dirs <- c("../../simul/fig_1/m_2000")
#angle = c(-0.2, 1.2)
angle = c(1,2,3) #c(-1.4, -1.3, -1.2 , -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4)
of        <- "fig1_modulopi"
modulo <- pi
#####################

df.all <- data.frame(NULL)

sims.dir <- list.dirs(sims.dirs)[2:length(list.dirs(sims.dirs))]
for (i in sims.dir) {
  df <- df.raster.map(i, modulo=modulo)
  pop <- str_split(i, "../../simul/fig_1/", n=2, simplify = TRUE)
  df[,8] <- sprintf("%s", pop[,2])
  df[,9] <- (df$ang_M-df$ang_S)
  df[10] <- modulo.all(df[,9], modulo=pi)
  df[11] <- df[,10]^2
  df.all <- rbind(df.all, df)
  
}

a <- -1.570796
b <- 0.03332000
modulo.all(a-b, modulo=pi)

a <- 1.570796
b <- 0.03332000
modulo.all(a-b, modulo=pi)

#Bug : for last directory, S angle is calculated as pi/1 and M angle as close to 0. Modulopi finds the biggest difference between the two, while it shouldn't be.
#Is it a problem of distance or of the calcul of S angle ?!! S is NOT pi/2 in those parameters file, it's the other way around. 
#Is there somewhere a mix up with A and B axis ? It doesn't happen when ellipses are drawn ; so I guesse it's in the measurement



#####################
sims.dirs <- c("../../simul/fig_1/w2")
#angle = c(-0.2, 1.2)
angle = c(1,2,3) #c(-1.4, -1.3, -1.2 , -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4)
of        <- "fig1_modulopi"
modulo <- pi
#####################

df.all <- data.frame(NULL)

sims.dir <- list.dirs(sims.dirs)[2:length(list.dirs(sims.dirs))]
for (i in sims.dir) {
  df <- df.raster.map(i, modulo=modulo)
  pop <- str_split(i, "../../simul/fig_1/", n=2, simplify = TRUE)
  df[,8] <- sprintf("%s", pop[,2])
  df[,9] <- (df$ang_M-df$ang_S)
  df[10] <- modulo.all(df[,9], modulo=pi)
  df[11] <- df[,10]^2
  df.all <- rbind(df.all, df)
  
}
