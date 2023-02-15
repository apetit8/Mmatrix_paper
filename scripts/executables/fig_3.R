source("scripts/functions_R/All_functions.R")
library(igraph)
#######################
sims.dirs <- c("simul/fig_2cd/0-grn")
#######################

df.topo <- df.data(sims.dirs, pattern = "simul/fig_2cd/", variable="netw", file_size=210000, w_of_6=TRUE, network=TRUE)

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



# plot(df.topo$a_b, df.topo$ang_M)
# plot( (df.topo$c_a+df.topo$a_c), df.topo$corrS)

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

pdf("figures/fig_3a.pdf", width=6, height=7)

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

pdf("figures/fig_3b.pdf", width=6, height=7)
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

pdf("figures/fig_3c.pdf", width=6, height=7)
plot.igraph(G, edge.width=(abs(E(G)$weight))^2*10, edge.color=ifelse(abs(E(G)$weight) < 0.4,"grey",ifelse(E(G)$weight > 0, "dodgerblue2","brown2")),
            layout =  layout_in_circle, edge.label=round(E(G)$weight,2), edge.label.color=ifelse(abs(E(G)$weight) > 0.4,"black","white"), edge.label.cex=1.5,
            vertex.label = c("d","b","a","c","e","f"), vertex.label.cex=2,
            vertex.size=25, vertex.label.color="black", edge.curved=0.2)
title(expression(paste("Correlation between \ngenes regulations and r(M)")), cex.main=1.5,line = -0.5)

dev.off()
