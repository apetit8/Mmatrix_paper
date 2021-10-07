source("../functions_R/network.R")
source("../functions_R/tools.R")
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(matrixStats)
#######################
sims.dirs <- list.dirs("../../simul/cercle/wag2/c0_0")
of        <- "NEWtopo_0_0"
where     <- "cercle"
sims.dir  <- sims.dirs[2:30]
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

M <- df[,1:12] #Dataframe

pp <- prcomp(M, scale = TRUE)


data.mean <- colMeans(M) # faster than apply(M, 2, mean)
data.sd <- apply(M, 2, sd)


# pdf(sprintf("../../figures/%s/%s_topoCP1.pdf", where, of), width=6, height=6) #PC1 inverse de celle utilisée pour tracer graph ..... Faire attention
par(mfrow = c(1, 3))

for (i in c(-5,0,5)) {
  
    myPCs <- matrix(c(i,0,0,0,0,0,0,0,0,0,0,0), ncol = 12)
    newWsc <- myPCs %*% t(pp$rotation)
    newW <- t( t(newWsc) * data.sd + data.mean)
    #unscale the data
    n <- (1+sqrt(1+4*ncol(M)))/2
    zeros <- n*(0:(n-1)) + 1:n
    neworder <- c(5,9,13,2,10,14,3,7,15,4,8,12,1,6,11,16) 
    #Messy, but used to draw matrix in the "right" way.
    #neworder <- c(which(!(1:(n*n)) %in% zeros), zeros) #
    newW <- cbind(newW, matrix(0, ncol=n, nrow=nrow(newW)))[,order(neworder)]
    listW <- as.matrix(rbind(c(newW[1:4]),c(newW[5:8]),c(newW[9:12]),c(newW[13:16])))
    #network drawing
    colnames(listW) <- NULL
    plot.network(listW, max=1.5)
    title(main=paste("PC1 = ", i), line=-2)
}

# dev.off()


# fviz_contrib(pp, choice = "var", axes = 1, top = 10)
# plot(pp, choix = "ind", autoLab = "yes")
# cp1 <- fviz_contrib(pp, choice = "var", axes = 1, top = 10)
# cp2 <- fviz_contrib(pp, choice = "var", axes = 2, top = 10)
# cp3 <- fviz_contrib(pp, choice = "var", axes = 3, top = 10)
# cp4 <- fviz_contrib(pp, choice = "var", axes = 4, top = 10)
# 
# grid.arrange(cp1, cp2, cp3, cp4, ncol=2)
# 
# 
# fviz_pca_ind(pp,
#              geom.ind = "point",
#              col.ind = df$M_angle, # colorer by groups
#              legend.title = "M angle", asp=1)+
#   scale_color_gradient2(low="darkblue",  mid="green", high="yellow", midpoint=mean(df$M_angle), space = "Lab")
# 
# 
