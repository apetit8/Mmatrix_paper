source("../functions_R/network.R")
source("../functions_R/tools.R")

#This script compares the quality of the M-S alignment between populations.

#######################
sims.dirs1 <- list.dirs("../../simul/fig_2/se_0_ad")
sims.dir1  <- sims.dirs1[2:length(sims.dirs1)]
#
sims.dirs2 <- list.dirs("../../simul/fig_2/se_0_da")
sims.dir2  <- sims.dirs2[2:length(sims.dirs2)]
#
of        <- "fig2"
modulo = pi
#######################

df.topo1 <- df.topo.raw(sims.dir1, network=FALSE)
df.topo2 <- df.topo.raw(sims.dir2, network=FALSE)

#T-test cannot be applied on data : angles does not follow a normal distribution. Moreover, the mean is meaningless.
#Kolmogorov-Smirnov test : compare distributions
ks.test(df.topo1$ang_M, df.topo2$ang_M)
#If the p-value <0.05, there is a statistical difference between the M angle distribution of the simulations.

diff1 <- (modulo.all(df.topo1$ang_M-df.topo1$ang_S, modulo=modulo))
diff2 <- (modulo.all(df.topo2$ang_M-df.topo2$ang_S, modulo=modulo))

#MSE comparison
mean(diff1^2)
hist(diff1^2)
mean(diff2^2)
hist(diff2^2)
#The pop with the lower mean has the best alignment.

#MSE comparison with Wilcoxon test (non parametric t-test as the distribution is not normal) :
wilcox.test((diff1^2),(diff2^2))
#A significant p-value indicate that there is a significant difference between the two alignment quality.

#Kolmogorov-Smirnov test : compare distributions (non parametric t-test as the distribution is not normal)
ks.test((diff1^2),(diff2^2))
#A significant p-value indicate that there is a significant difference between the two alignment quality.




# model1 <- lm(ang_M~A_B+A_C+A_D+B_A+B_C+B_D+C_A+C_B+C_D+D_A+D_B+D_C, data=df.topo1)
# model2 <- lm(ang_M~A_B+A_C+A_D+B_A+B_C+B_D+C_A+C_B+C_D+D_A+D_B+D_C, data=df.topo2)
# summary(model1)
# summary(model2)


