source("../functions_R/network.R")
source("../functions_R/tools.R")
#######################
sims.dirs <- c("../../simul/cercle/wag2/c-1_-1","../../simul/cercle/wag2/c0_0","../../simul/cercle/wag2/c1_1")
of        <- "PLS_ACP_topo2"
what      <- "angle" #"angle" or "ecc"
where     <- "cercle"
#######################

pdf(sprintf("../../figures/%s/%s_PLS1.pdf", where, of), width=5, height=5)
#layout(t(1:2))
#par(mfrow = c(1, 2))
for (i in sims.dirs) {
  
        sims.dir  <- list.dirs(i)[2:30]
      
      df.topo <- df.topo.raw(sims.dir)
      #topo regPLS
      
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
      df[] <- lapply( df, as.numeric)
      
      
      library(plsdepot)
      pls1 = plsreg1( df[,1:12], df[,13], comps = 2)
      summary(pls1)
      plot(pls1)
      #C'est les liaisons BA/AB qui influencent le + l'angle ?
      
      plot(df[,13], pls1$y.pred, type = "n", xlab="Original", ylab = "Predicted")
      title(sprintf("Comparison of responses of %s", i), cex.main = 0.9)
      abline(a = 0, b = 1, col = "gray85", lwd = 2)
      text(df[,13], pls1$y.pred, col = "#5592e3")
      
      #R2 for each of our components
      pls1$R2
      
}
dev.off()

