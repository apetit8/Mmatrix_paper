library(ggpubr)
library(ggplot2)
library(plyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/figure_tools.R")

#######################
sims.dirs <- list.dirs("../../simul/cercle/wag2/c1_1")
sims.dir  <- sims.dirs[2:30]
of        <- "small1_1"
what      <- "angle" #"angle" or "ecc"
where     <- "topo"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################


df.topo <- df.feat.and.topo(sims.dir, of, what, where, threshold, c)
#df.topo <- df.topo0.5_0.5


#ANOVA
df.aov <- as.data.frame( cbind(df.topo$what_M, df.topo$what_G, df.topo$matter, df.topo$what_S))

matter <- data.frame()
for (i in c(1:nrow(df.topo))) {
    if   (df.aov[i,3] != as.character("others") ){
      newt <- c(df.aov[i,])
    }
    else {
      newt <- c()
    }
    matter <- rbind(matter, newt)
}
# names(matter)[names(matter) == "V1"] <- "what_M"
# names(matter)[names(matter) == "V2"] <- "what_G"
# names(matter)[names(matter) == "V3"] <- "topo"
# names(matter)[names(matter) == "V4"] <- "what_S"
# df.aov <- as.data.frame( matter)

df.aov[,1] <- as.numeric(df.aov[,1])
df.aov[,2] <- as.numeric(df.aov[,2])
df.aov[,4] <- as.numeric(df.aov[,4])


# fit1 <- aov(what_M ~ topo, df.aov)
# 
# fit2 <- aov(what_G ~ topo, df.aov)

# sink((sprintf("../../figures/%s/anova/%s_aovangle_%s.txt", where, of, threshold)))
#     print("M what")
#     summary(fit1)
#     print("G what")
#     summary(fit2)
# sink()

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
nb.cols <- length(unique(df.aov$topo))
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)

if (what=="ecc") {lim=c(0,1)}
if (what=="angle") {lim=c(-pi/2,pi/2)}


pdf(sprintf("../../figures/%s/%s_topogrp_M.pdf", where, of), width=7, height=3)
ggplot(data=df.topo, aes(as.factor(round(x=what_M, digits = 1)) ,
                         y=frequency(what_M), Freq, fill = as.factor(matter)))+
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#02db77","#6ae8b1", "#8862ba", "#eb8698", "#E69F00", "#56B4E9", "#999999"))+
  theme(axis.text.x = element_text(angle = 90))+
  xlim(c(-1.5,1.5))+
  labs(title=sprintf("Topologies (%s) distribution in %s", length(unique(df.topo$topo)), of), x ="M angle", y = "Occurences /1160")
  
dev.off()


#BOXPLOT
pdfname7  <- print(sprintf("../../figures/%s/anova/%s_topos_AOV_%s.pdf", where, of, threshold))
pdf(pdfname7, width=15, height=4)
layout(1:1)

    #M on S
    ggplot(data=df.aov, aes(what_S, what_M))+
      coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
      geom_abline()+
      scale_color_manual(values = mycolors)+
      geom_point(pch=1, aes(color = topo))+
      labs(title=sprintf("Main topologies distribution in %s", of), x =what, y = "M what") 
    
    ggboxplot(df.aov, x = "topo", y ="what_M" , 
              color = "topo",
              ylab = "M_angle", xlab = "Topo")+
               scale_color_manual(values = mycolors)+
              theme(axis.text.x = element_text(angle = 90)) 
    par(mfrow = c(2, 2))
    plot(fit1)
    
    
    # #G on S
    # ggplot(data=df.aov, aes(what_S, what_G))+
    #   coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    #   geom_abline()+
    #   scale_color_manual(values = mycolors)+
    #   geom_point(pch=1, aes(color = topo))+
    #   labs(title=sprintf("Main topologies distribution in %s", of), x =what, y = "G what") 
    # 
    # ggboxplot(df.aov, x = "topo", y ="what_G" , 
    #           color = "topo",
    #           ylab = "G_angle", xlab = "Topo")+
    #           scale_color_manual(values = mycolors)+
    #           theme(axis.text.x = element_text(angle = 90)) 
    # par(mfrow = c(2, 2))
    # plot(fit2)
    
    #Topologies
    par(mfrow = c(3, 4))
    listW <- list()
    for (i in c(1:nrow(df.draw))) {
      lst <- (as.character(df.draw[i,3:18]))
      output <- matrix(as.numeric(unlist(lst)), ncol = 4, byrow = FALSE)
      listW <- list.append(listW,output)
    }
    
    for (i in unique(listW)) {
      draw.topos(i)
      title(main=paste( unlist(as.vector(t(i))), collapse=''), line=-2)
    }
    
dev.off()


#df.topo0.5_0.2_0.8 <-df.topo

