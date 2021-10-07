#multi for ggplot

source("../functions_R/network.R")
source("../functions_R/tools.R")


#######################
sims.dirs <- list.dirs("../../simul/cercle/wag2/c0_0")
sims.dir  <- sims.dirs[2:30]
of        <- "w2c0_0"
what      <- "angle" #"angle" or "ecc"
where     <- "wagner"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################

#features.extract.matrix(sims.dir)

df.multi <- df.feat.multi(sims.dir, mod="wagner")

# fit1 <- aov(ang_M ~ ang_S, df.multi)
# fit2 <- aov(ang_M ~ absdist_MS, df.multi)
# 
# sink((sprintf("../../figures/%s/%s_aovangle_%s.txt", where, of, threshold)))
#   print("M ang")
#   summary(fit1)
#   print("M Dist ang")
#   summary(fit2)
# sink()

pdf(sprintf("../../figures/%s/%s_angplot2_%s.pdf", where, of, threshold), width=8, height=5)
layout(t(1:7))
par(mfrow = c(3, 1))
#lim=c(0,1)
lim=c(-pi/2,pi/2) 

  ggplot(data=df.multi, aes())+
    coord_fixed(ratio = 1, xlim = lim, ylim = c(0,pi/2), expand = TRUE, clip = "on")+
    geom_point(pch=1, aes(ang_S, absdist_MS),  color="darkolivegreen3")+
    geom_smooth(aes(ang_S, absdist_MS), se=FALSE, color="darkolivegreen4")+
    labs(title=sprintf("M %s on S angle in %s", what, of), x ="S", y = "M")
  
  
  # ggboxplot(df.multi, x = "Value", y ="absdist_MS",
  #           ylab = "ang_M", xlab = "ang_S")+
  #   #scale_color_manual(values = mycolors)+
  #   labs(title=sprintf("M %s distance from S %s in %s", ang, ang, of), x ="S", y = "M")+
  #   theme(axis.text.x = element_text(angle = 90))
  # 
  # par(mfrow = c(2, 2))
  # plot(fit2)
  
  ggplot(data=df.multi, aes())+
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = TRUE, clip = "on")+
    geom_point(pch=1, aes(ang_S, ang_M),  color="darkolivegreen3")+
    # stat_summary(aes(ang_S, ang_M, color="orange"),fun=mean, geom="line")+
    geom_smooth(aes(ang_S, ang_M), se=FALSE, color="darkolivegreen4")+
    # geom_point(pch=1, aes(ang_S, ang_G, colour="G"))+
    # geom_smooth(aes(ang_S, ang_G, color='darkblue'),se=FALSE)+
    geom_abline()+
    labs(title=sprintf("M %s on S angle in %s", what, of), x ="S", y = "M")
  
  
  # ggboxplot(df.multi, x = "Value", y ="ang_M",
  #           ylab = "ang_M", xlab = "ang_S")+
  #   #scale_color_manual(values = mycolors)+
  #   labs(title=sprintf("M %s distribution on S %s in %s", ang, ang, of), x ="S", y = "M")+
  #   theme(axis.text.x = element_text(angle = 90)) 
  # par(mfrow = c(2, 2))
  # plot(fit1)
dev.off()
