source("../functions_R/network.R")
source("../functions_R/tools.R")
library(FactoMineR)
library(factoextra)
library(gridExtra)
#######################
sims.dirs <- list.dirs("../../simul/cercle/wag2/c-1_-1")
sims.dir  <- sims.dirs[2:30]
of        <- "TESTangle_plotacp_-1"
what      <- "angle" #"angle" or "ecc"
where     <- "cercle"
grp <- 5
scaling <- TRUE
#######################

df.topo <- df.topo.raw(sims.dir)


names(df.topo)[names(df.topo) == "V9"] <- "A-B"
names(df.topo)[names(df.topo) == "V10"] <- "A-C"
names(df.topo)[names(df.topo) == "V11"] <- "A-D"
names(df.topo)[names(df.topo) == "V12"] <- "B-A"
names(df.topo)[names(df.topo) == "V14"] <- "B-C"
names(df.topo)[names(df.topo) == "V15"] <- "B-D"
names(df.topo)[names(df.topo) == "V16"] <- "C-A"
names(df.topo)[names(df.topo) == "V17"] <- "C-B"
names(df.topo)[names(df.topo) == "V19"] <- "C-D"
names(df.topo)[names(df.topo) == "V20"] <- "D-A"
names(df.topo)[names(df.topo) == "V21"] <- "D-B"
names(df.topo)[names(df.topo) == "V22"] <- "D-C"

dfx <- df.topo[,c(9,10,11,12,14,15,16,17,19,20,21,22)]
dfy <- as.data.frame(df.topo$ang_M)
df <-as.data.frame(cbind(dfx,dfy))
names(df)[names(df) == "df.topo$ang_M"] <- "M_angle"
df[] <- lapply( df, as.numeric)
df.topo[] <- lapply( df.topo, as.numeric)


pdf(sprintf("../../figures/%s/%s_ACP.pdf", where, of), width=5, height=4)

# #Determine number of clusters
# wss <- (nrow(df[,1:12])-1)*sum(apply(df[,1:12],2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(df[,1:12],
#                                      centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares")
# 
# # K-Means Cluster Analysis
# fit <- kmeans(df[,1:12], grp) # grp cluster solution
# # get cluster means
# aggregate(df[,1:12],by=list(fit$cluster),FUN=mean)
# # append cluster assignment
# df <- data.frame(df, fit$cluster)
#   df$fit.cluster <- as.factor(df$fit.cluster)
#   df.topo <- data.frame(df.topo, df$fit.cluster)
# 
# # Ward Hierarchical Clustering
# d <- dist(df[,1:12], method = "euclidean") # distance matrix
# fit <- hclust(d, method="ward.D2")
# plot(fit) # display dendogram
# groups <- cutree(fit, k=grp) # cut tree into grp clusters
# # draw dendogram with red borders around the grp clusters
# rect.hclust(fit, k=grp, border="red")

  res.pca <- PCA(df[,1:12], scale.unit=scaling, graph = FALSE)
  eig.val <- res.pca$eig
  # fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))


  # cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
  # cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
  # cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
  # cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
  #
  # grid.arrange(cp1, cp2, cp3, cp4, ncol=2)
  #
#   plot(res.pca, choix = "var", autoLab = "yes")#correlation circle
#
# plot(res.pca, choix = "ind", autoLab = "yes") #plot PCA without groups

# # #PCA with groups
# fviz_pca_ind(res.pca,
#              geom.ind = "point", # Montre les points seulement (mais pas le "text")
#              col.ind = df$fit.cluster, # colorer by groups
#              #palette = c("#00AFBB", "#E7B800", "#FC4E07","#00CC00","#CC0066", "#808080", "#CC99FF", ),
#              addEllipses = TRUE, # Ellipses de concentration
#              legend.title = "Groups"
# )
# #PCA with angle
# fviz_pca_ind(res.pca,
#              geom.ind = "point",
#              col.ind = df$M_angle, # colorer by groups
#              legend.title = "M angle", asp=1)+
#              scale_color_gradient2(low="darkblue",  mid="green", high="yellow", midpoint=mean(df$M_angle), space = "Lab")
# 
# ggplot(data=df.topo, aes(ang_S, ang_M))+
#   coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
#   geom_abline()+
#   #scale_color_manual(values = mycolors)+
#   geom_point(pch=1, aes(color =as.factor(df.fit.cluster)))+
#   labs(title=sprintf("Main topologies distribution in %s", of), x =what, y = "M what")

# #Network drawing
# par(mfrow = c(2, 3))
# df.draw <- as.data.frame(cbind(as.factor(df.topo$df.fit.cluster), df.topo[,8:23]))
# for (j in unique(df.draw[,1])) {
#   print(j)
#   netw <- data.frame()
#   netw2 <- data.frame()
#   for (k in c(1:nrow(df.draw))) {
#     if   (as.character(df.draw[k,1]) == j ){
#       netw <- as.data.frame(c(df.draw[k,]))
#       }
#       else {      newt <- c()  }
#     netw2 <- as.data.frame(rbind(netw2, netw))
#     }
#   df.draw2 <- as.data.frame(netw2)
#   df.draw2[,2:17] <- lapply( df.draw2[,2:17], as.numeric)
#   df.net <- colMeans( df.draw2[,2:17])
# 
#   #listW <- list()
#   lst <- (as.character(df.net))
#   output <- matrix(as.numeric(unlist(lst)), ncol = 4, byrow = FALSE)
#   print(output)
#   # listW <- list.append(listW,output)
#   # draw.topos(output)
#   plot.network(output, max=0.5)
#   title(main=paste(j), line=-2)
# }

  #Plot angle with pc1 coordinates
  res.pca.values <- as.data.frame(res.pca$ind$coord)
  df.pc1 <- cbind(df.topo, res.pca.values)
  names(df.pc1)[names(df.pc1) == "Dim.1"] <- "pc1"
  names(df.pc1)[names(df.pc1) == "Dim.2"] <- "pc2"
  k <- which.max(df.pc1$pc1) 
  l <- which.min(df.pc1$pc1)
  
  library(viridis)
  ggplot(data=df.pc1, aes(ang_S, ang_M))+
    coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
    geom_abline()+
    scale_colour_gradientn(colours = c("darkblue", "green", "yellow"),limits = c(df.pc1[l,25], df.pc1[k,25]),na.value = "transparent")+
    #scale_colour_viridis(option="magma")+
    geom_point(pch=1, aes(color =pc1))+
    labs(title=sprintf("PC1 value distribution in %s", of), x ="Angle S", y = "Angle M")

  ggplot(data=df.pc1, aes(ang_S, ang_M))+
    coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
    geom_abline()+
    scale_colour_gradientn(colours = c("darkblue", "green", "yellow"),limits = c(df.pc1[l,25], df.pc1[k,25]),na.value = "transparent")+
    geom_point(pch=1, aes(color =pc2))+
    labs(title=sprintf("PC2 value distribution in %s", of), x ="Angle S", y = "Angle M")


dev.off()

# pdf(sprintf("../../figures/%s/%s_ACP_regression.pdf", where, of), width=5, height=4)
# for (i in df.pc1[,24:28]) {
#   plot(i,df.pc1$ang_M, main=cor(i,df.pc1$ang_M))
#   abline(lm(ang_M ~ i,df.pc1))
# }
# dev.off()

# # Networks draw "inside out" !! Correct before use
# pdf(sprintf("../../figures/%s/%s_topoCP1.pdf", where, of), width=6, height=6)
#     par(mfrow = c(1, 3))
#     #Draw networks
#     #Network CP max
#     k <- which.max(df.pc1$pc1)
#     output <- matrix(as.numeric(unlist(list(df.pc1[k,8:23]))), ncol = 4, byrow = FALSE)
#     output
#     plot.network(output, max=1)
#     title(main=paste("PC1 max = ",round(df.pc1[k,25], digits = 2)), line=-2)
#     #Network CP min
#     l <- which.min(df.pc1$pc1)
#     output <- matrix(as.numeric(unlist(list(df.pc1[l,8:23]))), ncol = 4, byrow = FALSE)
#     output
#     plot.network(output, max=1)
#     title(main=paste("PC1 min = ",round(df.pc1[l,25], digits = 2)), line=-2)
#     #Network CP mean
#     half <- mean(df.pc1[,25])
#     mean.n <- df.pc1[which(round(df.pc1$pc1, digits = 2) == round(half, digits = 2)), ]
#     output <- matrix(as.numeric(unlist(list(mean.n[1,8:23]))), ncol = 4, byrow = FALSE)
#     output
#     plot.network(output, max=1)
#     title(main=paste("PC1 moyenne= ", round(mean.n[1,25], digits = 2)), line=-2)
# 
# 
# dev.off()
