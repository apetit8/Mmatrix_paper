source("../functions_R/All_functions.R")
library(plsdepot)
#######################
sims.dirs <- c("../../simul/fig_2/2-full")
#######################

df.topo <- df.data(sims.dirs, pattern = "../../simul/fig_2/", variable="netw", file_size=220000, w_of_6=TRUE, network=TRUE)

dfx <- df.topo[,10:45]
dfy <- as.data.frame(df.topo$ang_M)
df <-as.data.frame(cbind(dfx,dfy))
df[] <- lapply( df, as.numeric)


df <- df.topo[,10:45]
df$alphaM <- df.topo$ang_M
df[] <- lapply( df, as.numeric)



pls1 = plsreg1( df[,c(2:7, 9:14, 16:21, 23:28, 30:35)], df[,37, drop = FALSE], comps = 2)
summary(pls1)
pls1$cor.xyt

# pdf("../../figures/PLS.pdf", width=7, height=7)
pdf("../../figures/fig2_PLS.pdf", width=6, height=6)
plot(pls1, what="variables", main = "", xlab = "Latent variable 1", ylab = "Latent variable 2")
dev.off()







###############################################################################
png(file="../../figures/PLS.png", width=400, height=400)
plot(pls1, what="variables", main = "")
dev.off()

g1 <- ggplot(df, aes(x=B_A, y=A_B, col=M_direction))+
  coord_fixed(ratio =1, ylim = c(-0.5,0.5), xlim = c(-0.5,0.5), expand = FALSE, clip = "on")+
  geom_point()+
  scale_color_gradientn(colours = c("yellow", "red", "blue", "yellow"), limits=c(-pi/2, pi/2), breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  ggtitle("B/ Regulation between A and B")+
  theme_bw()
