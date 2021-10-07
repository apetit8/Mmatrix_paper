library(ggplot2)
library(plyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/figure_tools.R")

#######################
sims.dirs <- list.dirs("../../simul/Wagner/topos/ecc1/l4_0.5_0.2-0.8_e1")
sims.dir  <- sims.dirs[2:8]
of        <- "l4_0.5_0.2-0.8"
what      <- "angle" #"angle" or "ecc"
where     <- "topo/ecc1"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################


df.topo <- df.feat.and.topo(sims.dir, of, what, where, threshold, c)

#palette
nb.cols <- length(unique(df.topo$complexity1))
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)

pdfname3  <- print(sprintf("../../figures/%s/%s_topos_compx_%s.pdf", where, of, threshold))
pdf(pdfname3, width=10, height=5)
layout((1:4))
    ggplot(as.data.frame(df.topo), aes(x=as.factor(df.topo$Value), y=frequency(df.topo$Value),  Freq, fill = as.factor(complexity1))) +
      scale_fill_manual(values = mycolors)+
      theme(axis.text.x = element_text(angle = 90)) +
      geom_bar(stat="identity") +
      labs(title=sprintf("Number of link in %s", of), x =what, y = "Effectif")
    
    ggplot(as.data.frame(df.topo), aes(x=as.factor(df.topo$Value), y=frequency(df.topo$Value),  Freq, fill = as.factor(complexity2))) +
      scale_fill_manual(values = mycolors)+
      theme(axis.text.x = element_text(angle = 90)) +
      geom_bar(stat="identity") +
      labs(title=sprintf("Number of link with non selected genes in %s", of), x =what, y = "Effectif")
    
    ggplot(as.data.frame(df.topo), aes(x=as.factor(df.topo$Value), y=frequency(df.topo$Value),  Freq, fill = as.factor(complexity3))) +
      scale_fill_manual(values = mycolors)+
      theme(axis.text.x = element_text(angle = 90)) +
      geom_bar(stat="identity") +
      labs(title=sprintf("Number of positive link in %s", of), x =what, y = "Effectif")
    
    ggplot(as.data.frame(df.topo), aes(x=as.factor(df.topo$Value), y=frequency(df.topo$Value),  Freq, fill = as.factor(complexity4))) +
      scale_fill_manual(values = mycolors)+
      theme(axis.text.x = element_text(angle = 90)) +
      geom_bar(stat="identity") +
      labs(title=sprintf("Number of negative in %s", of), x =what, y = "Effectif")
dev.off()

if (what=="ecc") {lim=c(0,1)}
if (what=="angle") {lim=c(-pi/2,pi/2)}

pdfname4  <- print(sprintf("../../figures/%s/%s_topo_compx-feat_%s.pdf", where, of, threshold))
pdf(pdfname4, width=10, height=5)
layout((1:4))
  #M on S
  ggplot(data=df.topo, aes(what_S, what_M))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity1)))+
    labs(title=sprintf("Number of link in %s", of), x =what, y = "M angle") 
  
  ggplot(data=df.topo, aes(what_S, what_M))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity2)))+
    labs(title=sprintf("Number of link with non selected genes in %s", of), x =what, y = "M angle") 
  
  ggplot(data=df.topo, aes(what_S, what_M))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity3)))+
    labs(title=sprintf("Number of positive link in %s", of), x =what, y = "M angle") 
  
  ggplot(data=df.topo, aes(what_S, what_M))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity4)))+
    labs(title=sprintf("Number of negative link in %s", of), x =what, y = "M angle") 
  
  #G on S
  ggplot(data=df.topo, aes(what_S, what_G))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity1)))+
    labs(title=sprintf("Number of link in %s", of), x =what, y = "G angle") 
  
  ggplot(data=df.topo, aes(what_S, what_G))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity2)))+
    labs(title=sprintf("Number of link with non selected genes in %s", of), x =what, y = "G angle")    
  
  ggplot(data=df.topo, aes(what_S, what_G))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity3)))+
    labs(title=sprintf("Number of positive link in %s", of), x =what, y = "G angle") 
  
  ggplot(data=df.topo, aes(what_S, what_G))+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity4)))+
    labs(title=sprintf("Number of negative link in %s", of), x =what, y = "G angle") 
dev.off()
