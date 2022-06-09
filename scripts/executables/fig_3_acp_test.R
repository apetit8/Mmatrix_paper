library(FactoMineR)
library(factoextra)
library(gridExtra)

df <- cbind(df.rnul$ang_S, df.rnul$ang_M, df.rnul$sq_dist, df.rnul$V12, df.rnul$V13, df.rnul$Fitness)

df <- as.data.frame(df)

names(df)[names(df) == "V1"] <- "ang_S"
names(df)[names(df) == "V2"] <- "ang_M"
names(df)[names(df) == "V3"] <- "d_alpha"
names(df)[names(df) == "V4"] <- "opt1"
names(df)[names(df) == "V5"] <- "opt2"
names(df)[names(df) == "V6"] <- "Fitness"



res.pca <- PCA(df, scale.unit=FALSE, graph = FALSE)
eig.val <- res.pca$eig

fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
grid.arrange(cp1, cp2, cp3, cp4, ncol=2)

plot(res.pca, choix = "var", autoLab = "yes")#correlation circle

fviz_pca_ind(res.pca,
             geom.ind = "point",
             col.ind = df$Fitness, # colorer by groups
             legend.title = "Fitness", asp=1)+
  scale_color_gradient2(low="darkblue",  mid="green", high="yellow", midpoint=mean(df$Fitness), space = "Lab")

