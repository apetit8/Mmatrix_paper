source("scripts/functions_R/All_functions.R")
library(plsdepot)
#######################
sims.dirs <- c("simul/fig_1de/3-grn")
#######################

df.topo <- df.data(sims.dirs, pattern = "simul/fig_1de/", variable="netw", file_size=210000, w_of_6=TRUE, network=TRUE)

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

# pdf("figures/PLS.pdf", width=7, height=7)
pdf("figures/fig2_PLS.pdf", width=6, height=6)
plot(pls1, what="variables", main = "", xlab = "Latent variable 1", ylab = "Latent variable 2")
dev.off()







###############################################################################
png(file="figures/PLS.png", width=400, height=400)
plot(pls1, what="variables", main = "")
dev.off()

g1 <- ggplot(df, aes(x=B_A, y=A_B, col=M_direction))+
  coord_fixed(ratio =1, ylim = c(-0.5,0.5), xlim = c(-0.5,0.5), expand = FALSE, clip = "on")+
  geom_point()+
  scale_color_gradientn(colours = c("yellow", "red", "blue", "yellow"), limits=c(-pi/2, pi/2), breaks=c(0, pi/4, pi/2, -pi/4, -pi/2),labels=c("0", "\u03c0/4", "\u03c0/2","-\u03c0/4", "-\u03c0/2"))+
  ggtitle("B/ Regulation between A and B")+
  theme_bw()


################################################################################
library(igraph)



W <- matrix(c(0,
              as.numeric(cor.test(df.topo$ang_M, df.topo$a_b)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$a_c)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$a_d)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$a_e)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$a_f)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$b_a)$estimate),0,
              as.numeric(cor.test(df.topo$ang_M, df.topo$b_c)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$b_d)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$b_e)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$b_f)$estimate), 
              as.numeric(cor.test(df.topo$ang_M, df.topo$c_a)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$c_b)$estimate), 0,
              as.numeric(cor.test(df.topo$ang_M, df.topo$c_d)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$c_e)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$c_f)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$d_a)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$d_b)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$d_c)$estimate),0,
              as.numeric(cor.test(df.topo$ang_M, df.topo$d_e)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$d_f)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$e_a)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$e_b)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$e_c)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$e_d)$estimate),0,
              as.numeric(cor.test(df.topo$ang_M, df.topo$e_f)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$f_a)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$f_b)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$f_c)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$f_d)$estimate),
              as.numeric(cor.test(df.topo$ang_M, df.topo$f_e)$estimate),0
              ), ncol = sqrt(length(df.topo[1,11:46])))



#W with the expected order for igraph
W1 <- W
W1[,1] <- W[,4]
W1[,3] <- W[,1]
W1[,4] <- W[,3]
W2 <- W1
W2[1,] <- W1[4,]
W2[3,] <- W1[1,]
W2[4,] <- W1[3,]
#


G <- as.directed(graph.adjacency(t(W2), weighted = T))
V(G)$color <- c("lightyellow", "lightpink", "lightpink", "lightyellow", "grey", "grey")
deg <- degree(G, mode = "all")
# G <- delete.edges(G, E(G)[ abs(weight) < 0.3 ])

pdf("figures/fig_graph_corr.pdf", width=6, height=7)

plot.igraph(G, edge.width=(abs(E(G)$weight))^2*10, edge.color=ifelse(abs(E(G)$weight) < 0.4,"grey",ifelse(E(G)$weight > 0, "dodgerblue2","brown2")),
            layout =  layout_in_circle, edge.label=round(E(G)$weight,2), edge.label.color=ifelse(abs(E(G)$weight) > 0.4,"black","white"), edge.label.cex=1.5,
            vertex.label = c("d","b","a","c","e","f"), vertex.label.cex=2,
            vertex.size=25, vertex.label.color="black", edge.curved=0.2)
title(expression(paste("Correlation between \ngenes regulations and ", alpha, "(M)")), cex.main=1.5,line = -0.5)

dev.off()


###ECC

W <- matrix(c(0,
              as.numeric(cor.test(df.topo$ecc_M, df.topo$a_b)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$a_c)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$a_d)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$a_e)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$a_f)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$b_a)$estimate),0,
              as.numeric(cor.test(df.topo$ecc_M, df.topo$b_c)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$b_d)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$b_e)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$b_f)$estimate), 
              as.numeric(cor.test(df.topo$ecc_M, df.topo$c_a)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$c_b)$estimate), 0,
              as.numeric(cor.test(df.topo$ecc_M, df.topo$c_d)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$c_e)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$c_f)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$d_a)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$d_b)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$d_c)$estimate),0,
              as.numeric(cor.test(df.topo$ecc_M, df.topo$d_e)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$d_f)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$e_a)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$e_b)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$e_c)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$e_d)$estimate),0,
              as.numeric(cor.test(df.topo$ecc_M, df.topo$e_f)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$f_a)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$f_b)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$f_c)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$f_d)$estimate),
              as.numeric(cor.test(df.topo$ecc_M, df.topo$f_e)$estimate),0
), ncol = sqrt(length(df.topo[1,11:46])))



#W with the expected order for igraph
W1 <- W
W1[,1] <- W[,4]
W1[,3] <- W[,1]
W1[,4] <- W[,3]
W2 <- W1
W2[1,] <- W1[4,]
W2[3,] <- W1[1,]
W2[4,] <- W1[3,]
#


G <- as.directed(graph.adjacency(t(W2), weighted = T))
V(G)$color <- c("lightyellow", "lightpink", "lightpink", "lightyellow", "grey", "grey")
deg <- degree(G, mode = "all")
# G <- delete.edges(G, E(G)[ abs(weight) < 0.1 ])

pdf("figures/fig_graph_corr_ecc.pdf", width=6, height=7)
plot.igraph(G, edge.width=(abs(E(G)$weight))^2*10, edge.color=ifelse(abs(E(G)$weight) < 0.1,"grey",ifelse(E(G)$weight > 0, "dodgerblue2","brown2")),
            layout =  layout_in_circle, edge.label=round(E(G)$weight,2), edge.label.color=ifelse(abs(E(G)$weight) > 0.1,"black","white"), edge.label.cex=1.5,
            vertex.label = c("d","b","a","c","e","f"), vertex.label.cex=2,
            vertex.size=25, vertex.label.color="black", edge.curved=0.2)
title(expression(paste("Correlation between \ngenes regulations and e(M)")), cex.main=1.5,line = -0.5)

dev.off()



###CORRS

W <- matrix(c(0,
              as.numeric(cor.test(df.topo$corrS, df.topo$a_b)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$a_c)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$a_d)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$a_e)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$a_f)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$b_a)$estimate),0,
              as.numeric(cor.test(df.topo$corrS, df.topo$b_c)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$b_d)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$b_e)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$b_f)$estimate), 
              as.numeric(cor.test(df.topo$corrS, df.topo$c_a)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$c_b)$estimate), 0,
              as.numeric(cor.test(df.topo$corrS, df.topo$c_d)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$c_e)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$c_f)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$d_a)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$d_b)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$d_c)$estimate),0,
              as.numeric(cor.test(df.topo$corrS, df.topo$d_e)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$d_f)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$e_a)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$e_b)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$e_c)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$e_d)$estimate),0,
              as.numeric(cor.test(df.topo$corrS, df.topo$e_f)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$f_a)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$f_b)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$f_c)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$f_d)$estimate),
              as.numeric(cor.test(df.topo$corrS, df.topo$f_e)$estimate),0
), ncol = sqrt(length(df.topo[1,11:46])))



#W with the expected order for igraph
W1 <- W
W1[,1] <- W[,4]
W1[,3] <- W[,1]
W1[,4] <- W[,3]
W2 <- W1
W2[1,] <- W1[4,]
W2[3,] <- W1[1,]
W2[4,] <- W1[3,]
#


G <- as.directed(graph.adjacency(t(W2), weighted = T))
V(G)$color <- c("lightyellow", "lightpink", "lightpink", "lightyellow", "grey", "grey")
deg <- degree(G, mode = "all")
# G <- delete.edges(G, E(G)[ abs(weight) < 0.1 ])

pdf("figures/fig_graph_corr_corrS.pdf", width=6, height=7)
plot.igraph(G, edge.width=(abs(E(G)$weight))^2*10, edge.color=ifelse(abs(E(G)$weight) < 0.4,"grey",ifelse(E(G)$weight > 0, "dodgerblue2","brown2")),
            layout =  layout_in_circle, edge.label=round(E(G)$weight,2), edge.label.color=ifelse(abs(E(G)$weight) > 0.4,"black","white"), edge.label.cex=1.5,
            vertex.label = c("d","b","a","c","e","f"), vertex.label.cex=2,
            vertex.size=25, vertex.label.color="black", edge.curved=0.2)
title(expression(paste("Correlation between \ngenes regulations and r(M)")), cex.main=1.5,line = -0.5)

dev.off()
