source("../functions_R/suppfig1_functions.R")
source("../functions_R/All_functions.R")
#####################
fig.width  <- 8
fig.height <- 4
n.points   <- 21
#####################

tW1 <- matrix(NA, nrow=6, ncol=6)
diag(tW1) <- 0

#A-B = 0
tW2 <- tW1
tW2[1,2] <- tW2[2,1] <- 0
#A-B = -0.5
tW3 <- tW1
tW3[1,2] <- tW3[2,1] <- -0.5
#A-B = 0.5
tW4 <- tW1
tW4[1,2] <- tW4[2,1] <- 0.5


oa1 <- optimW.alpharange(tW1, n.points=n.points, sigP=0.1)
pdf("../../figures/fig_supp1_A.pdf", width=fig.width, height=fig.height)
.plot_oa(oa1)
title("a-b regulations free to evolve", outer=TRUE, line=-3)
dev.off()

oa2 <- optimW.alpharange(tW2, n.points=n.points, sigP=0.1)
pdf("../../figures/fig_supp1_B.pdf", width=fig.width, height=fig.height)
  .plot_oa(oa2)
  title("a-b regulations fixed at 0", outer=TRUE, line=-3)
dev.off()

oa4 <- optimW.alpharange(tW4, n.points=n.points, sigP=0.1)
pdf("../../figures/fig_supp1_C.pdf", width=fig.width, height=fig.height)
.plot_oa(oa4)
title("a-b regulations fixed at +0.5", outer=TRUE, line=-3)
dev.off()

oa3 <- optimW.alpharange(tW3, n.points=n.points, sigP=0.1)
pdf("../../figures/fig_supp1_D.pdf",width=fig.width, height=fig.height)
  .plot_oa(oa3)
  title("a-b regulations fixed at -0.5", outer=TRUE, line=-3)
dev.off()


