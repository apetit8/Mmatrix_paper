source("../functions_R/All_functions.R")
library(plsdepot)
#######################
sims.dirs <- c("../../simul/fig_3/2-full")
of        <- "PLS_ACP_topo2"
what      <- "angle" #"angle" or "ecc"
where     <- "cercle"
#######################


df.topo <- df.data(sims.dirs, pattern = "../../simul/fig_2/", variable="netw", file_size=220000, w_of_6=TRUE, network=TRUE)

dfx <- df.topo[,10:45]
dfy <- as.data.frame(df.topo$ang_M)
df <-as.data.frame(cbind(dfx,dfy))
df[] <- lapply( df, as.numeric)


df <- df.topo[,10:45]
df$M_direction <- df.topo$ang_M
df[] <- lapply( df, as.numeric)



pls1 = plsreg1( df[,c(2:7, 9:14, 16:21, 23:28, 30:35)], df[,37, drop = FALSE], comps = 2)
summary(pls1)
pls1$cor.xyt

g1 <- ggplot(df, aes(x=B_A, y=A_B, col=M_direction))+
  coord_fixed(ratio =1, ylim = c(-0.5,0.5), xlim = c(-0.5,0.5), expand = FALSE, clip = "on")+
  geom_point()+
  scale_color_gradientn(colours = c("yellow", "red", "blue", "yellow"), limits=c(-pi/2, pi/2), breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  ggtitle("B/ Regulation between A and B")+
  theme_bw()

# pdf("../../figures/PLS.pdf", width=7, height=7)
pdf("../../figures/PLS.pdf", width=6, height=6)
  plot(pls1, what="variables", main = "PLS: Cercle of Correlations")
dev.off()


png(file="../../figures/PLS.png", width=400, height=400)
plot(pls1, what="variables", main = "PLS: Cercle of Correlations")
dev.off()


#######################
sims.dirs <- c("../../simul/fig_2/2-full")
#######################


df.topo <- df.data(sims.dirs, pattern = "../../simul/fig_2/", variable="netw", file_size=220000, w_of_6=TRUE, network=TRUE)

dfx <- df.topo[,10:45]
dfy <- as.data.frame(df.topo$corr)
df <-as.data.frame(cbind(dfx,dfy))
df[] <- lapply( df, as.numeric)


df <- df.topo[,10:45]
df$alpha_M <- df.topo$ang_M
df[] <- lapply( df, as.numeric)



pls1 = plsreg1( df[,c(2:7, 9:14, 16:21, 23:28, 30:35)], df[,37, drop = FALSE], comps = 2)
summary(pls1)
pls1$cor.xyt

g1 <- ggplot(df, aes(x=B_A, y=A_B, col=corr))+
  coord_fixed(ratio =1, ylim = c(-0.5,0.5), xlim = c(-0.5,0.5), expand = FALSE, clip = "on")+
  geom_point()+
  scale_color_gradientn(colours = c("yellow", "red", "blue", "yellow"), limits=c(-pi/2, pi/2), breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  ggtitle("B/ Regulation between A and B")+
  theme_bw()

# pdf("../../figures/PLS.pdf", width=7, height=7)
# pdf("../../figures/PLS.pdf", width=10, height=7)
# layout(t(1:2))
# par(mfrow = c(1, 2))
# for (i in sims.dirs) {


p <- plot(pls1, what="variables", main = "A/ PLS: Cercle of Correlations")
print(g1)
# 
# plot(df[,37], pls1$y.pred, type = "n", xlab="Original", ylab = "Predicted")
# title("Round S", cex.main = 0.9)
# abline(a = 0, b = 1, col = "gray85", lwd = 2)
# text(df[,37], pls1$y.pred, col = "#5592e3")
# #R2 for each of our components
# pls1$R2

# }
# dev.off()

