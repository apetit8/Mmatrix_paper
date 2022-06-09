#Script that atkes ALL info from output files for fig_3 format


df.raster.map <- function(sims.dir, modulo=pi, all.gen=FALSE){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  files     <- list.files(path = sims.dir, full.names=TRUE, pattern = "\\.txt$")
  #change here how to get to the right param file !
  files <- files[sapply(files, file.size) > 150000]
  for (i in files) { 
    #
    parname1<- str_split(i, ".txt$", n=2, simplify = TRUE)
    parname2<- str_split(parname1, "-R", n=2, simplify = TRUE)
    parname3<- str_split(parname2, "NEXTPAR-", n=2, simplify = TRUE)
    parname4 <- str_split(parname2, "/FITNESS_OPTIMUM-", n=2, simplify = TRUE)
    param.file <- sprintf("%s/%s", parname4[1,1], parname3[1,2]) #get the name of the second paramfile, with the right S
    #
    rp <- read.param(param.file)
    S.ref2 <- matrix.features(extract.S.matrix.fig3(param.file), n.genes=2)[["angle"]][1]
    
    mytopos <- lapply(i, function(ff) {
      tt <- read.table(ff, header=TRUE)
      mygens <-rev(if (all.gen==TRUE) tt[,"Gen"] else tt[nrow(tt),"Gen"])
      for (gen in mygens) {
          #M data
          gen <-gen
          phen.mean <- extract.P.mean(tt, gen=gen)
          M.mat <- extract.M.matrix(tt, gen=gen)
          M.ang_ppi <- matrix.features(M.mat,n.genes=2)[["angle"]][1]+pi #M matrix
          M.ang_mpi <- matrix.features(M.mat,n.genes=2)[["angle"]][1]-pi #M matrix
          M.feat1 <- matrix.features(M.mat,n.genes=2)[["eccentricity"]][1] #M matrix
          M.ang <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
          M.feat3 <- matrix.features(M.mat,n.genes=2)[["size"]][1] #M matrix
          fitness <- extract.fitness(tt, "MFit", gen=gen)
          Wmat <- extract.W.matrix(tt) #W data
          
          data.gen <- c(ff, gen, S.ref2, phen.mean[1],phen.mean[2],M.ang_ppi,M.ang_mpi, M.feat1, M.ang, M.feat3, fitness, Wmat )
          filedata <- rbind(filedata, data.gen)
      }
      return(as.list(filedata))
    })
    newt <- as.data.frame(rbindlist(mytopos, use.names=FALSE))
    setnames(newt, 1:27, c("data.dir","Gen","ang_S","P_mean_A","P_mean_B","ang_M0",
                           "ang_M2","ecc_M","ang_M","siz_M","Fitness",
                           "A_A","B_A","C_A","D_A","A_B","B_B","C_B","D_B","A_C",
                           "B_C","C_C","D_C","A_D","B_D","C_D","D_D"))
    

    # newt <- ldply(mytopos) # create a data.frame for mytopos
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  #Naming network cells (for 4 genes network)

  simul.df[,2:27] <- lapply( simul.df[,2:27], as.numeric)

  simul.df$sq_dist <- 1 - (((modulo.all(simul.df$ang_M-simul.df$ang_S, modulo = modulo))^2)/((pi^2)/12))
  return(simul.df)
}

