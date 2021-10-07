library(ggpubr)
library(ggplot2)
library(plyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(rlist)

source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/figure_tools.R")

#######################
sims.dirs <- list.dirs("../../simul/simu/l4_0.5_0.8")
sims.dir  <- sims.dirs[2:2]
of        <- "l4_0.5_0.2"
what      <- "angle" #"angle" or "ecc"
where     <- "topo"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################


df.topo <- df.feat.and.topo(sims.dir, of, what, where, threshold, c)




#ANOVA
df.dist <- as.data.frame( cbind(df.topo$absdist_MS, df.topo$absdist_GS, df.topo$matter, df.topo$what_S))
matter <- data.frame()
for (i in c(1:nrow(df.topo))) {
  if   (df.dist[i,3] != as.character("others") ){
    newt <- c(df.dist[i,])
  }
  else {
    newt <- c()
  }  
  matter <- rbind(matter, newt)
}
names(matter)[names(matter) == "V1"] <- "absdist_MS"
names(matter)[names(matter) == "V2"] <- "absdist_GS"
names(matter)[names(matter) == "V3"] <- "topo"
names(matter)[names(matter) == "V4"] <- "what_S"
df.dist <- as.data.frame( matter)

df.dist[,1] <- as.numeric(df.dist[,1])
df.dist[,2] <- as.numeric(df.dist[,2])
df.dist[,4] <- as.numeric(df.dist[,4])


fit1 <- aov(absdist_MS ~ topo, df.dist)
summary(fit1)


fit2 <- aov(absdist_GS ~ topo, df.dist)
summary(fit2)

sink((sprintf("../../figures/%s/anova/%s_AOVdist_%s.txt", where, of, threshold)))
print("M what")
summary(fit1)
print("G what")
summary(fit2)
sink()

####Draw
df.draw <- as.data.frame(cbind(df.topo$matter, df.topo$Value, df.topo[,5:20]))
matter <- data.frame()
for (i in c(1:nrow(df.draw))) {
  if   (df.draw[i,1] != as.character("others") ){
    newt <- c(df.draw[i,])  }
  else {
    newt <- c()  }  
  matter <- rbind(matter, newt)}
df.draw <- as.data.frame( matter)

#palette
nb.cols <- length(unique(df.dist$topo))
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)
if (what=="ecc") {lim=c(0,1)}
if (what=="angle") {lim=c(-pi/2,pi/2)}

#BOXPLOT
pdfname9  <- print(sprintf("../../figures/%s/anova/%s_topos_AOVdist_%s.pdf", where, of, threshold))
pdf(pdfname9, width=15, height=10)
layout(1:1)
    
    #M on S
    ggplot(data=df.dist, aes(what_S, absdist_MS))+
      coord_fixed(ratio = 1, xlim = lim, ylim = c(0,pi/2), expand = TRUE, clip = "on")+
      scale_color_manual(values = mycolors)+
      geom_point(pch=1, aes(color = topo))+
      labs(title=sprintf("Absolute distance between M and S in %s", of), x =what, y = "M absolute distance from S") 
    
    ggboxplot(df.dist, x = "topo", y ="absdist_MS" , 
              color = "topo",
              ylab = "M_angle", xlab = "Topo")+
      scale_color_manual(values = mycolors)+
      theme(axis.text.x = element_text(angle = 90)) 
    par(mfrow = c(2, 2))
    plot(fit1)
    
    
    #G on S
    ggplot(data=df.dist, aes(what_S, absdist_GS))+
      coord_fixed(ratio = 1, xlim = lim, ylim = c(0,pi/2), expand = TRUE, clip = "on")+
      scale_color_manual(values = mycolors)+
      geom_point(pch=1, aes(color = topo))+
      labs(title=sprintf("Absolute distance between G and S in %s", of), x =what, y = "G absolute distance from S") 
    
    ggboxplot(df.dist, x = "topo", y ="absdist_GS" , 
              color = "topo",
              ylab = "G_angle", xlab = "Topo")+
      scale_color_manual(values = mycolors)+
      theme(axis.text.x = element_text(angle = 90)) 
    par(mfrow = c(2, 2))
    plot(fit2)
    
    
    #Topologies
    par(mfrow = c(3, 4))
    listW <- list()
    for (i in c(1:nrow(df.draw))) {
      lst <- (as.character(df.draw[i,3:18]))
      output <- matrix(as.numeric(unlist(lst)), ncol = 4, byrow = TRUE)
      listW <- list.append(listW,output)
    }
    
    for (i in unique(listW)) {
      draw.topos(i)
      title(main=paste( unlist(as.vector(t(i))), collapse=''), line=-2)
    }
dev.off()





