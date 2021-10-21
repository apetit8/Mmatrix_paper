source("../functions_R/network.R")
source("../functions_R/tools.R")
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(matrixStats)
#######################
sims.dirs <- list.dirs("../../simul/cercle/wag2/c0_0")
of        <- "topo_0_1.5"
where     <- "cercle"
sims.dir  <- sims.dirs[2:30]
#######################

df.topo <- df.topo.raw(sims.dir)

names(df.topo)[names(df.topo) == "V9"] <- "A_B"
names(df.topo)[names(df.topo) == "V10"] <- "A_C"
names(df.topo)[names(df.topo) == "V11"] <- "A_D"
names(df.topo)[names(df.topo) == "V12"] <- "B_A"
names(df.topo)[names(df.topo) == "V14"] <- "B_C"
names(df.topo)[names(df.topo) == "V15"] <- "B_D"
names(df.topo)[names(df.topo) == "V16"] <- "C_A"
names(df.topo)[names(df.topo) == "V17"] <- "C_B"
names(df.topo)[names(df.topo) == "V19"] <- "C_D"
names(df.topo)[names(df.topo) == "V20"] <- "D_A"
names(df.topo)[names(df.topo) == "V21"] <- "D_B"
names(df.topo)[names(df.topo) == "V22"] <- "D_C"
df.topo[,9:23] <- lapply( df.topo[,9:23], as.numeric) 


#Linear regression
model <- lm(ang_M~ang_S, data=df.topo)

#Mean squared error :
data <-(data.frame(pred = df.topo$ang_S, actual = df.topo$ang_M))
data[] <- lapply( data, as.numeric) 
mean((data$actual - data$pred)^2)


#Linear regression with every W cell
model <- lm(ang_M~A_B+A_C+A_D+B_A+B_C+B_D+C_A+C_B+C_D+D_A+D_B+D_C, data=df.topo)
summary(model)

equation1=function(x){coef(model)[2]*x+coef(model)[1]}
equation2=function(x){coef(model)[2]*x+coef(model)[1]+coef(model)[3]}

ggplot(data=df.topo, aes(as.numeric(ang_S),as.numeric( ang_M)))+
  coord_fixed(ratio = 1, expand = TRUE, clip = "on")+
  geom_abline()+
  geom_point(aes(color =(A_B)))

ggplot(data=df.topo, aes(as.numeric(ang_M),as.numeric( A_B)))+
  geom_point()

ggplot(data=df.topo, aes(as.numeric(ang_M),as.numeric( B_A)))+
  geom_point()

ggplot(data=df.topo, aes(as.numeric(ang_M),as.numeric( C_D)))+
  geom_point()

ggplot(data=df.topo, aes(as.numeric(ang_M),as.numeric(A_C)))+
  geom_point()
