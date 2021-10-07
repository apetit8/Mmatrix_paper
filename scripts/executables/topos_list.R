source("../functions_R/network.R")
source("../functions_R/tools.R")
library(dplyr)
library(tidyr)
library(purrr)
#######################
sims.dirs <- list.dirs("../../simul/Wagner/topos/l4_0.5_0.35")
sims.dir  <- sims.dirs[2:30]
of        <- "l4_0.5_0.35"
where     <- "topo_csv"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################

j<-0

for (i in sims.dir) {
  files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
  j<- j+1
  
  
    mytopos <- lapply(files, function(ff) {
    #print(ff)
    tt <- read.table(ff, header=TRUE)
    Wmat <- extract.W.matrix(tt)
    cWmat <- cleanW(Wmat, a=c, epsilon=threshold) #use the threshold to group the topology who are considered the same
    untopo <-unique.topo(cWmat, groups=list(1:2,3:4))
    #print(paste(untopo, collapse="."))
  }
  )
  topos <- table(sapply(mytopos, function(topo) paste(topo, collapse=".")))
  topos.list <- data.frame(sort(topos, decreasing = TRUE))
  
  
  string <- stringr::str_remove(string = i, pattern = sims.dirs)
  csvname   <- print(sprintf("../../figures/topo_csv/%s%s_topos_list_%s.csv", of, string[1], threshold))
  write.table(data.frame(topos.list), file = csvname, col.names = TRUE)
  
}

print("Topos_list done !")