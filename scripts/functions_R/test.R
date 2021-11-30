mean.angle.2pi <- function(data) {
  atan2(mean(sin(data)), mean(cos(data)))
}

mean.angle.pi.byS <- function(data, ang_S) {
  df <- modulo.all(mean.angle.2pi(2*(data %% pi))/2, pi)
  if ((df - mean(ang_S, 7)/2) < -pi/2){ df = df + pi}
  if ((df - mean(ang_S, 7)/2) > pi/2){ df = df - pi}
  return(df)
}

modulo.all <- function(angle, modulo=pi) {
  angle <- angle %% modulo
  ifelse(angle > modulo/2, angle - modulo, angle)
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
      #M data
      gen <-tt[nrow(tt),"Gen"]
      M.mat <- extract.M.matrix(tt, gen=gen)
      M.feat1 <- matrix.features(M.mat,n.genes=2)[["eccentricity"]][1] #M matrix
      M.feat2 <- matrix.features(M.mat,n.genes=2)[["angle"]][1] #M matrix
      M.feat3 <- matrix.features(M.mat,n.genes=2)[["size"]][1] #M matrix
      mean.ma.df <- M.feat2
      if (network==TRUE){
        Wmat <- extract.W.matrix(tt) #W data
        new <- c(ff, S.ref1, S.ref2, S.ref3, M.feat1, M.feat2, M.feat3, mean.ma.df, Wmat )
      }
      if (network!=TRUE){
        new <- c(ff, S.ref1, S.ref2, S.ref3, M.feat1, M.feat2, M.feat3, mean.ma.df )
      }
      return(new)
      #create a list of data
    })
    newt <- ldply(mytopos) # create a data.frame for i
    newt[,8] <- mean.angle.pi.byS(data=as.numeric(newt[,8]), ang_S=as.numeric(newt[,3]))
    simul.df <- rbind(simul.df, newt) #paste the data of i in a data.frame with the others
  }
  simul.df[,2:8] <- lapply( simul.df[,2:8], as.numeric)
  
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

# Suffit d'ajouter + ang_S pour la moyenne ? Simple !
# Du coup non c'est pas ça mdr
#Il faudrait +pi dans certaines conditions.
#Conditions : if meanM - S/2 < -pi/2, alors meanM <- meanM + pi
#####################
sims.dirs <- c(
  "../../simul/fig_2/se_0","../../simul/fig_2/se_0_abba","../../simul/fig_2/se_0_o_abbacddc",
  "../../simul/fig_2/se_0_o_dacb"
)
of        <- "fig2"
modulo <- pi
#####################

df.all <- data.frame(NULL)
for (i in sims.dirs) {
  sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
  df <- df.topo.raw(sims.dir, network=FALSE)
  df[,9] <- modulo.all(df$ang_M-df$ang_S, modulo = modulo)
  sq_rho <- round(mse.ms(df[,9]), 4)
  pop <- str_split(i, "../../simul/fig_2/", n=2, simplify = TRUE)
  df[,10] <- sprintf("%s", pop[,2])
  df.all <- rbind(df.all, df)
}
names(df.all)[names(df.all) == "V9"] <- "ang_diff"
names(df.all)[names(df.all) == "V10"] <- "netw"

p2 <- ggplot(data=df.all, aes(ang_S, ang_M,col=netw))+
  coord_fixed(ratio = 1)+
  geom_abline()+
  labs(title = ("A/ M angle distribution and mean for different S angle"), y="M Angle", x="S Angle")+
  # geom_smooth(aes(fill=netw), alpha = 0, method = "AngleSmooth", span=0.9, level = 0.95)+
  geom_point(aes(y=mean_ang_M, fill=netw), show.legend = FALSE)+
  stat_summary( aes(ang_S, ang_M,col=netw), geom = "linerange", position = "identity", fun.data = st_dev_angle, size=3, alpha=0.2, show.legend = FALSE)
p2 <- fracAx(p=p2, y=TRUE, "pi")
p2 <- p2 + facet_wrap(vars(netw), ncol=2) + theme(strip.background = element_blank(), strip.text = element_blank())
p2
