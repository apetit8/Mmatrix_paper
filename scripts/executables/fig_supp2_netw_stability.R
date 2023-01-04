#Fig supp
source("scripts/functions_R/All_functions.R")
#####################
sims.dirs <- c(list.dirs("simul/fig_3", recursive = FALSE), "simul/fig_1de/3-grn", list.dirs("simul/fig_4", recursive = FALSE))
#####################

df.stab <- df.data(sims.dirs, pattern = "simul/fig_3/", variable="netw", file_size=15000, w_of_6=TRUE, network=TRUE)

df.stability <- data.frame()
for (i in 1:nrow(df.stab)) {
  W <- (matrix(as.numeric(df.stab[i,11:(ncol(df.stab)-1)]), ncol = 6) )
  
  st <- data.frame(
    data.dir   = df.stab[i,1], 
    pop = df.stab[i,ncol(df.stab)],
    Gen= df.stab[i,2],
    n.step = stablefrom(W))
  
  df.stability <- rbind(df.stability, st)
}

hist(df.stability$n.step, breaks=6)


cairo_pdf("figures/figsupp_stability.pdf", width=4, height=5)
hist(df.stability$n.step, main = NULL, freq = FALSE, breaks =seq(1, 24, 1),
     xlab = "Timesteps before gene expression stability", ylab = "Frequency", xlim = c(3,24))

dev.off()

