mean.angle.2pi <- function(angles) {
  atan2(mean(sin(angles)), mean(cos(angles)))
}

df.topo.raw<- function(sims.dir, w_of_4=TRUE, network=TRUE){
  simul.df <- data.frame()

  #collect all the data
  for (i in sims.dir) {
    files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
    param.file = list.files(path=i, pattern="param.par", full.names=TRUE)
    stopifnot(length(files) >0, length(param.file) == 1)

    #read param.par to associate each simu with the S "what"
    rp <- read.param(param.file)
    mut.rate <- 2*rp$GENET_MUTRATES
    if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    S.ref1 <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["eccentricity"]][1]
    S.ref2 <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["angle"]][1]
    S.ref3 <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["size"]][1]

    mytopos <- lapply(files, function(ff) {
      #print(ff)
      tt <- read.table(ff, header=TRUE)
      if (network==TRUE){
        Wmat <- extract.W.matrix(tt)
      }
      #"what" of G and M
      gen <-tt[nrow(tt),"Gen"]
      P.mat <- extract.P.matrix(tt, gen=gen)
      M.mat <- extract.M.matrix(tt, gen=gen)
      M.feat1 <- matrix.features(M.mat,n.genes=2)[["eccentricity"]][1] #M matrix
      M.feat2 <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
      M.feat3 <- matrix.features(M.mat,n.genes=2)[["size"]][1] #M matrix
      mean.ma.df <- M.feat2
      if (network==TRUE){
        new <- c(ff, S.ref1, S.ref2, S.ref3, M.feat1, M.feat2, M.feat3, mean.ma.df, Wmat )
      }
      if (network!=TRUE){
        new <- c(ff, S.ref1, S.ref2, S.ref3, M.feat1, M.feat2, M.feat3, mean.ma.df )
      }
      return(new)
      #create a list of data
    })
    newt <- ldply(mytopos) # create a data.frame for i
    newt[,8] <- mean.angle.2pi(as.numeric(newt[,8]))
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  simul.df[,2:8] <- lapply( simul.df[,2:8], as.numeric)

  # if (mean_m==TRUE){
  #   s_ang <- unique(round(simul.df[,3],3))
  #   meandf <- data.frame()
  #   for (k in s_ang){
  #     dfdf <- filter(simul.df[,3:6], round(simul.df[,3],3) == k)
  #     browser()
  #     dfdf[,4] <- mean.angle.2pi(dfdf[,3])
  #     meandf %>% add_row(.dfdf, dfdf[,4])
  #     # meandf <- rbind(unname(meandf), data.frame(unname(dfdf[,4])))
  # 
  #       }
  #   simul.df[,8] <- meandf[,1]
  # }
  names(simul.df)[names(simul.df) == "V1"] <- "data.dir"
  names(simul.df)[names(simul.df) == "V2"] <- "ecc_S"
  names(simul.df)[names(simul.df) == "V3"] <- "ang_S"
  names(simul.df)[names(simul.df) == "V4"] <- "siz_S"
  names(simul.df)[names(simul.df) == "V5"] <- "ecc_M"
  names(simul.df)[names(simul.df) == "V6"] <- "ang_M"
  names(simul.df)[names(simul.df) == "V7"] <- "siz_M"
  names(simul.df)[names(simul.df) == "V8"] <- "mean_ang_M"

  if (w_of_4==TRUE & network==TRUE) {
    #Naming network cells (for 4 genes network)
    names(simul.df)[names(simul.df) == "V10"] <- "B_A"
    names(simul.df)[names(simul.df) == "V11"] <- "C_A"
    names(simul.df)[names(simul.df) == "V12"] <- "D_A"
    names(simul.df)[names(simul.df) == "V13"] <- "A_B"
    names(simul.df)[names(simul.df) == "V15"] <- "C_B"
    names(simul.df)[names(simul.df) == "V16"] <- "D_B"
    names(simul.df)[names(simul.df) == "V17"] <- "A_C"
    names(simul.df)[names(simul.df) == "V18"] <- "B_C"
    names(simul.df)[names(simul.df) == "V20"] <- "D_C"
    names(simul.df)[names(simul.df) == "V21"] <- "A_D"
    names(simul.df)[names(simul.df) == "V22"] <- "B_D"
    names(simul.df)[names(simul.df) == "V23"] <- "C_D"
    simul.df[,9:24] <- lapply( simul.df[,9:24], as.numeric)
  }

  return(simul.df)
}
