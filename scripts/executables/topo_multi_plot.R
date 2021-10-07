library(ggplot2)
library(plyr)
library(tidyr)
library(tidyverse)

source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/figure_tools.R")


#######################
sims.dirs <- list.dirs("../../simul/Wagner/topos/l4_0.5_0.35")
sims.dir  <- sims.dirs[2:30]
of        <- "l4_0.5_0.35"
what      <- "angle" #"angle" or "ecc"
where     <- "topo"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################

df.topo <- df.feat.and.topo(sims.dir, of, what, where, threshold, c)

ggplot(as.data.frame(df.topo), aes(x=as.factor(df.topo$Value), y=frequency(df.topo$Value),  Freq, fill = as.factor(complexity1))) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_bar(stat="identity") +
  labs(title=of, x =what, y = "Effectif")

ggplot(as.data.frame(df.topo), aes(x=as.factor(df.topo$Value), y=frequency(df.topo$Value),  Freq, fill = what_G)) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_bar(stat="identity") +
  labs(title=of, x =what, y = "Effectif")

ggplot(data = df.topo, aes(x = what_S, y = what_G, color = complexity1)) +
  #geom_line() + 
  #geom_point() +
  geom_abline(intercept = 37, slope = -5)+
  ggtitle(eq)
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title=of, x =what, y = "Effectif")
  
  ggplot(data=df.topo, aes(what_S, what_M, col=complexity1))+
    geom_point(aes(col=complexity1))+
    stat_summary(aes(col=complexity1),fun=mean, geom="line")
  
  df2 <- cbind(df.topo[,2:4],df.topo$Value)
  df2 <- melt(df2 ,  id.vars = 'what_S', variable.name = 'series')
  ggplot(df2, aes(what_S,value)) + 
    geom_point(aes(colour = series)) +
    #geom_smooth(aes(fill='series'), span=0.75)
  
  ggplot(data=df.topo, aes(what_S, what_G, what_M))+
    geom_point(aes(col=complexity1))+
    stat_summary(aes(col=complexity1),fun=mean, geom="line")
  
  
  
  
