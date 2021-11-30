



df.raster.map <- function(sims.dir, modulo=pi, all.gen=FALSE){
  simul.df <- data.frame()
  filedata <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    files     <- list.files(path = i, full.names=TRUE, pattern = "\\.txt$")
    param.file.all = list.files(path=i, pattern="\\.par$", full.names=TRUE)
    param.file <- as.character(param.file.all[1])
    stopifnot(length(files) >0)
    
    #read param.par to associate each simu with S properties
    rp <- read.param(param.file)
    mut.rate <- 2*rp$GENET_MUTRATES
    if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    S.ref2 <- matrix.features(extract.S.matrix(param.file), n.genes=2)[["angle"]][1]
    

      mytopos <- lapply(files, function(ff) {
        #print(ff)
        tt <- read.table(ff, header=TRUE) 
        mygens <-rev(if (all.gen==TRUE) tt[,"Gen"] else tt[nrow(tt),"Gen"])
        for (gen in mygens) {
            phen.mean <- extract.P.mean(tt, gen=gen)
            M.mat <- extract.M.matrix(tt, gen=gen)
            M.feat2 <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
            data.gen <- c(ff, gen, S.ref2, M.feat2, phen.mean[1],phen.mean[2] )
            filedata <- rbind(filedata, data.gen)
            #create a list of data
        }
        return(as.list(filedata))
      })
      newt <- rbindlist(mytopos, use.names=FALSE)
      setnames(newt, 1:6, c("data.dir", "Gen","ang_S","ang_M","P_mean_A","P_mean_B"))
      # newt <- ldply(mytopos) # create a data.frame for mytopos
      simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others

  }
  simul.df[,2:6] <- lapply( simul.df[,2:6], as.numeric)
  sq_dist <- (modulo.all(simul.df$ang_M-simul.df$ang_S, modulo = modulo))^2
  simul.df <- as.data.frame(cbind(simul.df, sq_dist))
  return(simul.df)
  
}


# mygens <-rev(if (all.gen) simuls.mean[,"Gen"] else simuls.mean[nrow(simuls.mean),"Gen"])
#   for (gen in mygens) {
#     mytopos <- lapply(files, function(ff) {
#       #print(ff)
#       tt <- read.table(ff, header=TRUE) 
#       gen <-tt[nrow(tt),"Gen"]
#       phen.mean <- extract.P.mean(tt, gen=gen)
#       M.mat <- extract.M.matrix(tt, gen=gen)
#       M.feat2 <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
#       new <- c(ff, S.ref2, M.feat2, phen.mean[1],phen.mean[2] )
#       return(new)
#       #create a list of data
#     })
#     newt <- ldply(mytopos) # create a data.frame for i
#     simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
#   }
# }