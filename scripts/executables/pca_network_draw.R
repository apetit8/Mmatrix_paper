source("../functions_R/network.R")
source("../functions_R/tools.R")
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(matrixStats)
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
                "../../simul/wagner/4g_fp_ab","../../simul/wagner/4g_no_ab","../../simul/wagner/4g_chn_ab",
                "../../simul/wagner/4g_chp_ab","../../simul/wagner/4g_fhn_ab","../../simul/wagner/4g_fhp_ab"
                )
of        <- "controlled_W"
where     <- "wagner"
sims.dir  <- sims.dirs[2:30]
#######################

pdf(sprintf("../../figures/%s/%s_ACP.pdf", where, of), width=6, height=6)
for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:30]
    
    df.topo <- df.topo.raw(sims.dir)
    dfx <- df.topo[,c(9,10,11,12,14,15,16,17,19,20,21,22)]
    dfy <- as.data.frame(df.topo$ang_M)
    df <-as.data.frame(cbind(dfx,dfy))
    names(df)[names(df) == "df.topo$ang_M"] <- "M_angle"
    df[] <- lapply( df, as.numeric)
    df.topo[] <- lapply( df.topo, as.numeric)
    
    #Naming network cells (for 4 genes network)
    names(df.topo)[names(df.topo) == "V9"] <- "B-A"
    names(df.topo)[names(df.topo) == "V10"] <- "C-A"
    names(df.topo)[names(df.topo) == "V11"] <- "D-A"
    names(df.topo)[names(df.topo) == "V12"] <- "A-B"
    names(df.topo)[names(df.topo) == "V14"] <- "C-B"
    names(df.topo)[names(df.topo) == "V15"] <- "D-B"
    names(df.topo)[names(df.topo) == "V16"] <- "A-C"
    names(df.topo)[names(df.topo) == "V17"] <- "B-C"
    names(df.topo)[names(df.topo) == "V19"] <- "D-C"
    names(df.topo)[names(df.topo) == "V20"] <- "A-D"
    names(df.topo)[names(df.topo) == "V21"] <- "B-D"
    names(df.topo)[names(df.topo) == "V22"] <- "C-D"
    
    #extraction of the data for the PCA
    dfx <- df.topo[,c(9,10,11,12,14,15,16,17,19,20,21,22)]
    dfy <- as.data.frame(df.topo$ang_M)
    df <-as.data.frame(cbind(dfx,dfy))
    names(df)[names(df) == "df.topo$ang_M"] <- "M_angle"
    df[] <- lapply( df, as.numeric)  
    M <- df[,1:12] #Dataframe
    
    #PCA
    pp <- prcomp(M, scale = TRUE)
    ppp <- fviz_screeplot(pp, addlabels = TRUE, ylim = c(0, 50), title = i)
    print(ppp)

      #PC component
      cp1 <- fviz_contrib(pp, choice = "var", axes = 1, top = 10)
      cp2 <- fviz_contrib(pp, choice = "var", axes = 2, top = 10)
      cp3 <- fviz_contrib(pp, choice = "var", axes = 3, top = 10)
      cp4 <- fviz_contrib(pp, choice = "var", axes = 4, top = 10)
      grid.arrange(cp1, cp2, cp3, cp4, ncol=2)
      
      #Plot individuals in PCA, colored by M angle
      plpc <- fviz_pca_ind(pp,
                   geom.ind = "point",
                   col.ind = df$M_angle, # colors by M angle
                   legend.title = "M angle", asp=1)+
        scale_color_gradient2(low="yellow",  mid="green", high="darkblue", midpoint=mean(df$M_angle), space = "Lab")
      print(plpc)
      
      data.mean <- colMeans(M) # faster than apply(M, 2, mean)
      data.sd <- apply(M, 2, sd)
      
        #Plot angle with pc coordinates
        res.ind <- get_pca_ind(pp)
        res.pca.values <- as.data.frame(res.ind$coord)
        df.pc1 <- cbind(df.topo, res.pca.values)
        names(df.pc1)[names(df.pc1) == "Dim.1"] <- "pc1"
        names(df.pc1)[names(df.pc1) == "Dim.2"] <- "pc2"
        k <- max(df.pc1$pc1, na.rm=TRUE)
        l <- min(df.pc1$pc1, na.rm=TRUE)
        #pc1
        pl1 <- ggplot(data=df.pc1, aes(ang_S, ang_M))+
          coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
          geom_abline()+
          scale_colour_gradientn(colours = c("darkblue", "green", "yellow"),limits = c(l, k),na.value = "transparent")+
          geom_point(pch=1, aes(color =pc1))+
          labs(title=sprintf("PC1 value distribution in %s", of), x ="Angle S", y = "Angle M")
        print(pl1)
        #pc2
        pl2 <- ggplot(data=df.pc1, aes(ang_S, ang_M))+
          coord_fixed(ratio = 1, xlim = c(-pi/2,pi/2), ylim = c(-pi/2,pi/2), expand = TRUE, clip = "on")+
          geom_abline()+
          scale_colour_gradientn(colours = c("darkblue", "green", "yellow"),limits = c(l, k),na.value = "transparent")+
          geom_point(pch=1, aes(color =pc2))+
          labs(title=sprintf("PC2 value distribution in %s", of), x ="Angle S", y = "Angle M")
        print(pl2)
        
        #Network associated with PC1 values
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
        
}   
dev.off()

