library(ggplot2)
library(plyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

source("../functions_R/network.R")
source("../functions_R/tools.R")
source("../functions_R/figure_tools.R")


#######################
sims.dirs <- list.dirs("../../simul/Wagner/topos/l4_0.5_0.5")
sims.dir  <- sims.dirs[2:30]
of        <- "l4_0.5_0.5"
what      <- "angle" #"angle" or "ecc"
where     <- "topo/ecc0.12"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################

df.topo <- df.feat.and.topo(sims.dir, of, what, where, threshold, c)

#palette
nb.cols <- length(unique(df.topo$complexity1))
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)


pdfname  <- print(sprintf("../../figures/%s/%s_topos-n-feat_%s.pdf", where, of, threshold))
pdf(pdfname, width=10, height=5)
layout(t(1:2))

  #M on S
if (what == "ecc") {
  ggplot(data=df.topo, aes(what_S, what_M))+
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1), expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity1)))+
    labs(title=of, x =what, y = "M eccentricity")  
}
if (what == "angle") { 
  ggplot(data=df.topo, aes(what_S, what_M))+
    coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
    geom_abline()+
    scale_color_manual(values = mycolors)+
    geom_point(aes(color = as.factor(complexity1)))+
    labs(title=of, x =what, y = "M angle") 
}

#   #G on S
# #M on S
# if (what == "ecc") {
#   ggplot(data=df.topo, aes(what_S, what_G))+
#     coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1), expand = TRUE, clip = "on")+
#     geom_abline()+
#     scale_color_manual(values = mycolors)+
#     geom_point(aes(color = as.factor(complexity1)))+
#     labs(title=of, x =what, y = "G angle") 
# }
# if (what == "angle") { 
#   ggplot(data=df.topo, aes(what_S, what_G))+
#     coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-1.5,1.5), expand = TRUE, clip = "on")+
#     geom_abline()+
#     scale_color_manual(values = mycolors)+
#     geom_point(aes(color = as.factor(complexity1)))+
#     labs(title=of, x =what, y = "G angle") 
# }

  # #distance G on S
  # ggplot(data=df.topo, aes(color = as.factor(complexity1)))+
  #   scale_color_manual(values = mycolors)+
  #   geom_point(aes(x = Value, y = dist_GS))+
  #   #geom_point(pch=2, aes(x = Value, y = dist_MS))+
  # labs(title=of, x =what, y = "distance from S, G")  
  
  #distance M on S
  ggplot(data=df.topo, aes(color = as.factor(complexity1)))+
    scale_color_manual(values = mycolors)+
    geom_point(aes(x = Value, y = dist_MS))+
    #geom_point(pch=2, aes(x = Value, y = dist_MS))+
    labs(title=of, x =what, y = "distance from S, M")  
  
  # #absolute distance G on S
  # ggplot(data=df.topo, aes(color = as.factor(complexity1)))+
  #   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
  #   scale_color_manual(values = mycolors)+
  #   geom_point(aes(x = Value, y = absdist_GS))+
  #   #geom_point(pch=2, aes(x = Value, y = absdist_MS))+
  #   labs(title=of, x =what, y = "absolute distance from S, G") 

  #absolute distance M on S
  ggplot(data=df.topo, aes(color = as.factor(complexity1)))+
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
    scale_color_manual(values = mycolors)+
    geom_point(aes(x = Value, y = absdist_MS))+
    #geom_point(pch=2, aes(x = Value, y = absdist_MS))+
    labs(title=of, x =what, y = "absolute distance from S, M") 
dev.off()

