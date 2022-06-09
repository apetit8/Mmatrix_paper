source("../functions_R/All_functions.R")
#####################
sims.dirs <- c(
  "../../simul/fig_3_v2/no_ab"
)
modulo <- pi
#####################

dfv2.rnul <- df.opt.map2(sims.dirs, modulo=modulo, gen = FALSE)

dir <- str_split(dfv2.rnul$data.dir, ".par-R", n=3, simplify = TRUE)
dfv2.rnul$dir2 <- dir[,1]

dfv2_angevolv <- data.frame()
j <-1
for (u in unique(as.list(dir[,1]))){
  dfv2f <- subset(dfv2.rnul, dir2 == u)
  dfv2f[12] <- mean(dfv2f$sq_dist)
  dfv2_angevolv <- rbind(dfv2_angevolv, dfv2f[12] )
  j <-j+1
}
dfv2.rnul$d_alpha <- dfv2_angevolv[,1] #d_alpha

dir <- str_split(dfv2.rnul$dir2, ".FILE", n=3, simplify = TRUE)
dir <- str_split(dir[,1], "OPTIMUM-", n=3, simplify = TRUE)
dfv2.rnul$opti <- as.numeric(dir[,2]) #Optimum

# write.csv(dfv2.rnul, "../data/fig3-no_ab_fullgen_fullnetw.csv")

dfv2.rnul5000 <- subset(dfv2.rnul, Gen==5000)

dfv2.rnul <- dfv2.rnul5000



gn1 <- ggplot(data=dfv2.rnul)+
  stat_summary(aes(opti, B_A), fun=mean, geom="line", color="red", show.legend = FALSE)+
  stat_summary(aes(opti, A_B), fun=mean, geom="line", color="red", show.legend = FALSE)+
  #
  stat_summary(aes(opti, B_A), fun.data = mean_sdl, geom="errorbar", color="red", show.legend = FALSE)+
  stat_summary(aes(opti, A_B), fun.data = mean_sdl, geom="errorbar", color="red", show.legend = FALSE)+
  labs(y="Regulation", x="Optimum", title="Direct regulation between selected genes")  

gn2 <- ggplot(data=dfv2.rnul)+
  stat_summary(aes(opti, C_A), fun=mean, geom="line", color="darkblue", show.legend = FALSE)+
  stat_summary(aes(opti, D_A), fun=mean, geom="line", color="darkblue", show.legend = FALSE)+
  stat_summary(aes(opti, C_B), fun=mean, geom="line", color="darkblue", show.legend = FALSE)+
  stat_summary(aes(opti, D_B), fun=mean, geom="line", color="darkblue", show.legend = FALSE)+
  #
  stat_summary(aes(opti, C_A), fun.data = mean_sdl, geom="errorbar", color="darkblue", show.legend = FALSE)+
  stat_summary(aes(opti, D_A), fun.data = mean_sdl, geom="errorbar", color="darkblue", show.legend = FALSE)+
  stat_summary(aes(opti, C_B), fun.data = mean_sdl, geom="errorbar", color="darkblue", show.legend = FALSE)+
  stat_summary(aes(opti, D_B), fun.data = mean_sdl, geom="errorbar", color="darkblue", show.legend = FALSE)+
  labs(y="Regulation", x="Optimum", title="Direct regulation of non-selected genes to selected genes")  

gn3 <- ggplot(data=dfv2.rnul)+
  stat_summary(aes(opti, A_C), fun=mean, geom="line", color="forestgreen", show.legend = FALSE)+
  stat_summary(aes(opti, B_C), fun=mean, geom="line", color="forestgreen", show.legend = FALSE)+
  stat_summary(aes(opti, A_D), fun=mean, geom="line", color="forestgreen", show.legend = FALSE)+
  stat_summary(aes(opti, B_D), fun=mean, geom="line", color="forestgreen", show.legend = FALSE)+
  #
  stat_summary(aes(opti, A_C), fun.data = mean_sdl, geom="errorbar", color="forestgreen", show.legend = FALSE)+
  stat_summary(aes(opti, B_C), fun.data = mean_sdl, geom="errorbar", color="forestgreen", show.legend = FALSE)+
  stat_summary(aes(opti, A_D), fun.data = mean_sdl, geom="errorbar", color="forestgreen", show.legend = FALSE)+
  stat_summary(aes(opti, B_D), fun.data = mean_sdl, geom="errorbar", color="forestgreen", show.legend = FALSE)+
  labs(y="Regulation", x="Optimum", title="Direct regulation of selected genes to non-selected genes")  

gn4 <- ggplot(data=dfv2.rnul)+
  stat_summary(aes(opti, D_C), fun=mean, geom="line", color="purple", show.legend = FALSE)+
  stat_summary(aes(opti, C_D), fun=mean, geom="line", color="purple", show.legend = FALSE)+
  #
  stat_summary(aes(opti, D_C), fun.data = mean_sdl, geom="errorbar", color="purple", show.legend = FALSE)+
  stat_summary(aes(opti, C_D), fun.data = mean_sdl, geom="errorbar", color="purple", show.legend = FALSE)+
  labs(y="Regulation", x="Optimum", title="Direct regulation between non-selected genes")  


cairo_pdf("../../figures/fig_opt_reg.pdf", width=10, height=7)
grid.arrange(
  gn1, gn2, gn3, gn4,
  ncol = 2,
  nrow = 2,
  widths = c(1, 1),
  clip = FALSE
)
dev.off()
