source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs <- list.dirs("../../simul/fig_light_constraint", recursive = FALSE)
modulo <- pi
#####################


df.topo.wide.m<- function(sims.dir, w_of_4=FALSE, w_of_6=FALSE, network=FALSE, file_size=100000){
  simul.df <- data.frame()
  #collect all the data
  for (i in sims.dir) {
    print(i)
    files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
    files <- files[sapply(files, file.size) > file_size]
    # paramS <- str_split(i, "0.5/", n=2, simplify = TRUE)[2]
    paramSdir <- str_split(i, "simu-", n=2, simplify = TRUE)
    paramS <- str_split(paramSdir[2], ".par-nextpar", n=2, simplify = TRUE)[1]
    param.file = list.files(path=i, pattern="param.par", full.names=TRUE)
    param.fileS = paste0(paramSdir[1],paramS)
    # stopifnot(length(files) >0, length(param.file) == 1)
    #read param.par to associate each simu with the S "what"
    rp <- read.param(param.fileS)
    # mut.rate <- 2*rp$GENET_MUTRATES
    # if (rp$GENET_MUTTYPE=="locus") mut.rate <- mut.rate*rp$GENET_NBLOC
    S.ref2 <- matrix.features(M=extract.S.matrix.fig3(param.fileS), n.genes=2)[["angle"]][1] #S matrix
    
    mytopos <- lapply(files, function(ff) {
      print(ff)
      tt <- read.table(ff, header=TRUE)
      #M data
      gen <-tt[nrow(tt),"Gen"]
      M.mat <- extract.M.matrix(tt, gen=gen)
      M.ang_ppi <- matrix.features(M.mat,n.genes=2)[["angle"]][1]+pi #M matrix
      M.ang_mpi <- matrix.features(M.mat,n.genes=2)[["angle"]][1]-pi #M matrix
      M.feat1 <- matrix.features(M.mat,n.genes=2)[["eccentricity"]][1] #M matrix
      M.ang <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
      M.feat3 <- matrix.features(M.mat,n.genes=2)[["size"]][1] #M matrix
      fitness <- extract.fitness(tt, "MFit", gen=gen)
      opt1 <- extract.fitness(tt, "FitOpt1", gen=gen)[1]
      opt2 <- extract.fitness(tt, "FitOpt2", gen=gen)[1]
      if (network==TRUE){
        Wmat <- extract.W.matrix(tt) #W data
        new <- c(ff, M.ang_ppi, S.ref2, M.ang_mpi, M.feat1, M.ang, M.feat3, fitness, opt1, opt2, Wmat )
      }
      if (network!=TRUE){
        new <- c(ff, M.ang_ppi, S.ref2, M.ang_mpi, M.feat1, M.ang, M.feat3, fitness, opt1, opt2 )
      }
      return(new)
      #create a list of data
    })
    newt <- ldply(mytopos) # create a data.frame for i
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  
  if (network!=TRUE){
    simul.df <- setnames(simul.df[,1:10], c("data.dir","ang_M_ppi","ang_S","ang_M_mpi","ecc_M","ang_M","siz_M","fitness","opti_A","opti_B"))
    simul.df[,2:10] <- lapply( simul.df[,2:10], as.numeric)
  }
  if (w_of_6==TRUE & network==TRUE) {
    simul.df <- setnames(simul.df[,1:46], c("data.dir","ang_M_ppi","ang_S","ang_M_mpi","ecc_M","ang_M","siz_M","fitness","opti_A","opti_B",
                                            "A_A","B_A","C_A","D_A","E_A","F_A","A_B","B_B","C_B","D_B","E_B","F_B","A_C","B_C","C_C","D_C","E_C","F_C",
                                            "A_D","B_D","C_D","D_D","E_D","F_D","A_E","B_E","C_E","D_E","E_E","F_E","A_F","B_F","C_F","D_F","E_F","F_F"))
    simul.df[,2:46] <- lapply( simul.df[,2:46], as.numeric)
  }
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
    simul.df[,9:25] <- lapply( simul.df[,9:25], as.numeric)
  }
  return(simul.df)
}


#Data
df.fig2 <- df.data(sims.dirs, pattern = "../../simul/fig_light_constraint/", variable="netw", file_size=15000, w_of_6=TRUE, network=FALSE)

#With eccentricity
pfig2 <- ggplot(data=df.fig2, aes(ang_S, ang_M))+
  coord_fixed(ratio = 1, xlim = c(-1.5,1.5), ylim = c(-pi/2-0.2,pi/2+0.2), expand = TRUE, clip = "on")+
  geom_abline(colour="#666666")+
  geom_abline(intercept=pi, colour="#666666")+
  geom_abline(intercept=-pi, colour="#666666")+
  geom_point(aes(y=ang_M, col=1-ecc_M), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_mpi, col=1-ecc_M), alpha=0.16, show.legend = FALSE)+
  geom_point(aes(y=ang_M_ppi, col=1-ecc_M), alpha=0.16, show.legend = TRUE)+
  labs(y="M direction", x="S direction", fill = expression("\u03BE\u03B1"))+
  scale_color_viridis_c(option = "plasma")+
  labs(col = "M Eccentricity")+
  scale_x_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  scale_y_continuous(breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),
                     labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))
pfig2 <- pfig2 + facet_wrap(pop ~., ncol=3) + theme(strip.background = element_blank())+ theme_bw() #, strip.text = element_blank()
pfig2