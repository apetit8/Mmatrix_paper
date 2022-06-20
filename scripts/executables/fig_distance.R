source("../functions_R/All_functions.R")
library(png)
library(igraph)
#####################
sims.dirs <-  c(list.dirs("../../simul/fig_4", recursive = FALSE),list.dirs("../../simul/fig_3", recursive = FALSE))
modulo <- pi
#####################
# sims.dirs <-  c("../../simul/fig_2/2-full")
sims.dirs <-  list.dirs("../../simul/fig_4", recursive = FALSE)
df.fig4 <- df.data(sims.dirs, pattern = "../../simul/fig_4/", variable="netw", file_size=15000, w_of_6=TRUE, network=TRUE)

for (i in 1:nrow(df.fig4)){
  W <- t(matrix(as.numeric(df.fig4[i,10:45]), ncol = 6))
  W[,1:6]<-ifelse(abs(W[,1:6])>0.2,1,W[,1:6])
  W[,1:6]<-ifelse(abs(W[,1:6])<0.2,0,W[,1:6])
  df.fig4[i,]$siz_M <- mean(distances(graph_from_adjacency_matrix(as.table(W)))[1,2],
                            distances(graph_from_adjacency_matrix(as.table(W)))[1,2])
  # df.fig4[i,]$siz_M <-abs(W[1,2])+abs(W[2,1])
}

ggplot(df.fig4, aes(x=as.factor(siz_M), y=corr))+
  geom_violin(aes(fill=as.factor(siz_M)), show.legend = FALSE)+ theme_bw()+
  labs(y="Mutational co-expression", x="Mean connection distance with network")

ggplot(df.fig4, aes(x=as.factor(pop), y=corr, col=as.factor(pop)))+
  geom_point() #YES

ggplot(df.fig4, aes(x=as.factor(pop), y=corr))+
  geom_boxplot(aes(fill=as.factor(pop)), show.legend = FALSE)+ theme_bw()+
  labs(y="Mutational co-expression", x="Mean connection distance with network")

ggplot(df.fig4, aes(x=as.factor(siz_M), y=corr))+
  geom_boxplot(aes(fill=as.factor(siz_M)), show.legend = FALSE)+ theme_bw()+
  labs(y="Mutational co-expression", x="Mean connection distance with network")




library(ggstatsplot)
library(tidyverse)
plt <- ggbetweenstats(
  data = df.fig4,
  x = pop,
  y = corr,
  plot.type = "box",
  ggtheme = ggplot2::theme_bw()+theme(plot.margin = unit(c(0, 0, 3.3, 0), "cm")),
  pairwise.comparisons=FALSE,
  bf.message=FALSE,
  results.subtitle=FALSE,
  xlab="Network", ylab="Mutational coexpression",
  title = "Mutational co-expression evolved with different distance between A and B"
)
plt


png.netw1 = readPNG('../../templates/fig_4/networks_fig4.png')

cairo_pdf("../../figures/fig_distances.pdf", width=7, height=5.8)
grid.arrange(
  plt,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
grid.raster(png.netw1, x=0.54, y=0.13, width=0.9)
dev.off()



# W <- t(matrix(as.numeric(df.fig2[155,10:45]), ncol = 6))
# plotnet(W=W, gene.names=c("A","B","C","D","E","F"), thresh=0, max=1)
# 
# W[,1:6]<-ifelse(abs(W[,1:6])>0.1,1,W[,1:6])
# W[,1:6]<-ifelse(abs(W[,1:6])<0.1,0,W[,1:6])
# W
# 
# graph <- graph_from_adjacency_matrix(as.table(W))
# plot(graph)
# 
# mean(distances(graph)[1,2], distances(graph)[1,2])
# 

