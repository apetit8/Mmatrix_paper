source("../functions_R/network.R")
source("../functions_R/tools.R")


sims.dirs2 <- list.dirs("../../simul/fig_3/se_1")
sims.dir2  <- sims.dirs2[2:30]

df.topo2 <- df.topo.raw(sims.dir2)

df <- as.data.frame(cbind(df.topo2$ang_M,df.topo2$ang_S))
df[,1] <- df[,1]*-1

# ggplot(data = df.topo2, aes(ang_S,ang_M)) +
#   geom_point()
# 
# ggplot(data = df, aes(V2,V1)) +
#   geom_point()

ff <-((df[,1] - df[,2]))
dd <- modulo.mse(ff,modulo=pi/2)

mean((modulo.mse(ff,modulo=pi/2))^2)

#MSE comparison
mean(modulopi.diff(df[,1] - df[,2])^2)
hist(modulopi.diff(df[,1] - df[,2])^2)
mean(modulopi.diff(df.topo2$ang_M - df.topo2$ang_S)^2)
hist(modulopi.diff(df.topo2$ang_M - df.topo2$ang_S)^2)


mean(angle.diff(df[,1],df[,2]))
hist(angle.diff(df[,1],df[,2]))
mean(angle.diff(df.topo2$ang_M,df.topo2$ang_S))
hist(angle.diff(df.topo2$ang_M,df.topo2$ang_S))


#MSE comparison with Wilcoxon test (non parametric t-test as the distribution is not normal) :
wilcox.test((angle.diff(df[,1],df[,2])^2),(angle.diff(df.topo2$ang_M,df.topo2$ang_S)^2))
#A significant p-value indicate that there is a significant difference between the two alignment quality.

#appears more appropriate
#Kolmogorov-Smirnov test : compare distributions (non parametric t-test as the distribution is not normal)
ks.test((angle.diff(df[,1],df[,2])^2),(angle.diff(df.topo2$ang_M,df.topo2$ang_S)^2))
#A significant p-value indicate that there is a significant difference between the two alignment quality.


