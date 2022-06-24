source("../functions_R/plotnet.R")
source("../functions_R/All_functions.R")
library(png)
#####################
sims.dirs.fig2c <-list.dirs(
  "../../simul/fig_2_bis", recursive = FALSE
)

# sims.dirs.fig2c <-c("../../simul/fig_2_bis/0" , "../../simul/fig_2_bis/1", "../../simul/fig_2_bis/2", "../../simul/fig_2_bis/3","../../simul/fig_2_bis/4"                    )
modulo <- pi
#####################

df.datas <- data.frame(NULL)
for (i in sims.dirs.fig2c) {
  sims.dir <- list.dirs(i)[2:length(list.dirs(i))]
  df <- df.topo.wide.m(sims.dir, w_of_6=TRUE, network=TRUE)
  df.datas <- rbind(df.datas, df)
}



W <- matrix(as.numeric(df.datas[1,9:44]), ncol = 6)
plotnet(W=W, gene.names=c("A","B","C","D","E","F"), thresh=0, max=1)


W <- matrix(as.numeric(df.optis[300,11:46]), ncol = 6)
plotnet(W=W, gene.names=c("A","B","C","D","E","F"), thresh=0, max=1.5)
