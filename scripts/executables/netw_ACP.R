source("../functions_R/network.R")
source("../functions_R/tools.R")
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(viridis)
#######################
sims.dirs <- c(     
                # "../../simul/wagner/4g/c-1_-1","../../simul/wagner/4g/c1_1","../../simul/wagner/4g/c0_0",
                # "../../simul/wagner/8g/c-1_-1","../../simul/wagner/8g/c1_1","../../simul/wagner/8g/c0_0",
                # "../../simul/wagner/12g/c-1_-1","../../simul/wagner/12g/c1_1","../../simul/wagner/12g/c0_0",
                # "../../simul/wagner/16g/c-1_-1","../../simul/wagner/16g/c1_1","../../simul/wagner/16g/c0_0",
                # "../../simul/wagner/4g_pop500/c0_0","../../simul/wagner/4g_pop500/c1_1","../../simul/wagner/4g_pop500/c-1_-1",
                # "../../simul/wagner/4g_pop10000/c0_0","../../simul/wagner/4g_pop10000/c1_1","../../simul/wagner/4g_pop10000/c-1_-1"
                "../../simul/wagner/4g_cn_ab","../../simul/wagner/4g_cp_ab","../../simul/wagner/4g_fn_ab",
                "../../simul/wagner/4g_fp_ab","../../simul/wagner/4g_no_ab"
              )
of        <- "controlled_W"
what      <- "angle" #"angle" or "ecc"
where     <- "wagner"
grp <- 5
scaling <- TRUE
#######################

pdf(sprintf("../../figures/%s/%s_ACP.pdf", where, of), width=5, height=4)

for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:30]
      df.topo <- df.topo.raw(sims.dir)

      dfx <- df.topo[,c(9,10,11,12,14,15,16,17,19,20,21,22)]
      dfy <- as.data.frame(df.topo$ang_M)
      df <-as.data.frame(cbind(dfx,dfy))
      names(df)[names(df) == "df.topo$ang_M"] <- "M_angle"
      df[] <- lapply( df, as.numeric)
      df.topo[] <- lapply( df.topo, as.numeric)
      
      
        res.pca <- PCA(df[,1:12], scale.unit=scaling, graph = FALSE)
        eig.val <- res.pca$eig
        fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), title = i)

      
        cp1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
        cp2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
        cp3 <- fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
        cp4 <- fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
      
        grid.arrange(cp1, cp2, cp3, cp4, ncol=2)
      
      #   plot(res.pca, choix = "var", autoLab = "yes")#correlation circle
      #
      # plot(res.pca, choix = "ind", autoLab = "yes") #plot PCA without groups
      
      #PCA with angle
      fviz_pca_ind(res.pca,
                   geom.ind = "point",
                   col.ind = df$M_angle, # colorer by groups
                   legend.title = "M angle", asp=1)+
                   scale_color_gradient2(low="darkblue",  mid="green", high="yellow", midpoint=mean(df$M_angle), space = "Lab")
      

        #Plot angle with pc1 coordinates
        res.pca.values <- as.data.frame(res.pca$ind$coord)
        df.pc1 <- cbind(df.topo, res.pca.values)
        names(df.pc1)[names(df.pc1) == "Dim.1"] <- "pc1"
        names(df.pc1)[names(df.pc1) == "Dim.2"] <- "pc2"
        k <- max(df.pc1$pc1, na.rm=TRUE)
        l <- min(df.pc1$pc1, na.rm=TRUE)
        
        library(viridis)
        ggplot(data=df.pc1, aes(ang_S, ang_M))+
          coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
          geom_abline()+
          scale_colour_gradientn(colours = c("darkblue", "green", "yellow"),limits = c(l, k),na.value = "transparent")+
          #scale_colour_viridis(option="magma")+
          geom_point(pch=1, aes(color =pc1))+
          labs(title=sprintf("PC1 value distribution in %s", of), x ="Angle S", y = "Angle M")
      
        ggplot(data=df.pc1, aes(ang_S, ang_M))+
          coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
          geom_abline()+
          scale_colour_gradientn(colours = c("darkblue", "green", "yellow"),limits = c(l, k),na.value = "transparent")+
          geom_point(pch=1, aes(color =pc2))+
          labs(title=sprintf("PC2 value distribution in %s", of), x ="Angle S", y = "Angle M")
}

dev.off()

# pdf(sprintf("../../figures/%s/%s_ACP_regression.pdf", where, of), width=5, height=4)
# for (i in df.pc1[,24:28]) {
#   plot(i,df.pc1$ang_M, main=cor(i,df.pc1$ang_M))
#   abline(lm(ang_M ~ i,df.pc1))
# }
# dev.off()
