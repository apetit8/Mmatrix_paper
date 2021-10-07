source("../functions_R/network.R")
source("../functions_R/tools.R")
library(ggplot2)
library(plyr)
library(tidyr)
library(rlist)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)

#######################
sims.dirs <- list.dirs("../../simul/Wagner/topo_ecc1/l4_0.5_0.8")
sims.dir  <- sims.dirs[2:8]
of        <- "l4_0.5_0.8"
what      <- "angle" #"angle" or "ecc"
where     <- "topo/ecc1"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################


df.topo <- df.feat.and.topo(sims.dir, of, what, where, threshold, c)
#df.topo <- df.topo0.5_0.8

#palette
nb.cols <- length(unique(df.topo$matter))
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)
if (what=="ecc") {lim=c(0,1)}
if (what=="angle") {lim=c(-pi/2,pi/2)}   

pdfname6  <- print(sprintf("../../figures/%s/anova/%s_topo_draw_%s.pdf", where, of, threshold))
pdf(pdfname6, width=8, height=5)
layout(t(1:7))
  par(mfrow = c(2, 1))

    ggplot(data=df.topo, aes(what_S, what_M))+
      coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
      geom_abline()+
      scale_color_manual(values = mycolors)+
      geom_point(pch=1, aes(color = matter))+
      labs(title=sprintf("Main topologies distribution in %s", of), x =what, y = "M what")
    
    par(mfrow = c(1, 7))
      sims.dir  <- sims.dirs[4]
      draw.topo.angle.from.dir(sims.dir, threshold)
      sims.dir  <- sims.dirs[3]
      draw.topo.angle.from.dir(sims.dir, threshold)
      sims.dir  <- sims.dirs[2]
      draw.topo.angle.from.dir(sims.dir, threshold)
      sims.dir  <- sims.dirs[7]
      draw.topo.angle.from.dir(sims.dir, threshold)
      sims.dir  <- sims.dirs[5]
      draw.topo.angle.from.dir(sims.dir, threshold)
      sims.dir  <- sims.dirs[7]
      draw.topo.angle.from.dir(sims.dir, threshold)
      sims.dir  <- sims.dirs[8]
      draw.topo.angle.from.dir(sims.dir, threshold)

dev.off()

#dftopo0.2_8 <- df.topo


df.draw <- as.data.frame(cbind(df.topo$matter, df.topo$Value, df.topo[,5:20]))

####
layout(t(1:5))
matter <- data.frame()
for (i in c(1:nrow(df.draw))) {
  if   (df.draw[i,1] != as.character("others") ){
    newt <- c(df.draw[i,])
  }
  else {
    newt <- c()
  }
  matter <- rbind(matter, newt)
}
df.draw <- as.data.frame( matter)

#draw each topologies
listW <- list()
for (i in c(1:nrow(df.draw))) {
  lst <- (as.character(df.draw[i,3:18]))
  output <- matrix(as.numeric(unlist(lst)), ncol = 4)
  listW <- list.append(listW,output)
}

for (i in unique(listW)) {
  draw.topos(i)
  title(main=paste( unlist(i), collapse=''), line=-2)
}
