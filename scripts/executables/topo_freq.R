library(ggplot2)
library(plyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/figure_tools.R")

#######################
sims.dirs <- list.dirs("../../simul/cercle/wag2/c0_0")
sims.dir  <- sims.dirs[2:30]
of        <- "wag2small0_0"
what      <- "angle" #"angle" or "ecc"
where     <- "topo"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################


df.topo <- df.feat.and.topo(sims.dir, of, what, where, threshold, c)
#df.topo0.35 <- df.topo

#add frequency column
freq <- data.frame()
for (i in df.topo$topo) {
  occ <- df.topo[which(df.topo$topo == i), ]
  freq <- rbind(freq, nrow(occ))
}
df.topo$Frequency <- freq


matter <- data.frame()
for (i in c(1:nrow(df.topo))) {
  if (df.topo[i,32] >= nrow(df.topo)*0.015) {
    newt <- c(df.topo[i,21])
  }
  else {
    newt <- c("others")
  }  
  matter <- rbind(unname(matter), unname(newt))
}
df.topo$matter <- matter


#palette
nb.cols <- length(unique(df.topo$matter))
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)

if (what=="ecc") {lim=c(0,1)}
if (what=="angle") {lim=c(-pi/2,pi/2)}


pdfname5  <- print(sprintf("../../figures/%s/%s_topo_freq_%s.pdf", where, of, threshold))
pdf(pdfname5, width=10, height=5)
layout((1:4))
#M on S
  ggplot(data=df.topo, aes(what_S, what_M))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    #scale_fill_manual(values = mycolors)+
    geom_point(aes(color = Frequency))+
    labs(title=sprintf("Topology frequency in %s", of), x =what, y = "M what") 
 
  ggplot(data=df.topo, aes(what_S, what_M))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(pch=1, aes(color = matter))+
    labs(title=sprintf("Main topologies distribution in %s", of), x =what, y = "M what") 

# #G on S
#   ggplot(data=df.topo, aes(what_S, what_G))+
#     coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
#     geom_abline()+
#     #scale_color_manual(values = mycolors)+
#     geom_point(aes(color = Frequency))+
#     labs(title=sprintf("Topology frequency in %s", of), x =what, y = "M what") 
#   
#   ggplot(data=df.topo, aes(what_S, what_G))+
#     coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
#     geom_abline()+
#     scale_color_manual(values = mycolors)+
#     geom_point(pch=1, aes(color = matter))+
#     labs(title=sprintf("Main topologies distribution in %s", of), x =what, y = "M what") 

dev.off()






  
  