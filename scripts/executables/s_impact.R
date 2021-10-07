source("../functions_R/network.R")
source("../functions_R/tools.R")


#######################
sims.dirs1 <- list.dirs("../../simul/cercle/wag2/c0_0")
sims.dir1  <- sims.dirs1[2:30]
#
sims.dirs2 <- list.dirs("../../simul/cercle/wag2/c-1_-1")
sims.dir2  <- sims.dirs2[2:30]
#
sims.dirs3 <- list.dirs("../../simul/cercle/wag2/c1_1")
sims.dir3  <- sims.dirs3[2:30]
#
of        <- "wecc2_opt"
where     <- "wagner"
#######################

df.1 <- df.feat.multi(sims.dir1, mod="optimum_0")
df.2 <- df.feat.multi(sims.dir2, mod="optimum_-1")
df <- rbind(df.1, df.2)
df.3 <- df.feat.multi(sims.dir3, mod="optimum_1")
df <- rbind(df, df.3)

pdf(sprintf("../../figures/%s/%s_s_impact.pdf", where, of), width=8, height=5)
    layout(t(1:1))
    
    ggplot(data=df, aes(ang_S, absdist_MS, col=model))+
      geom_point(aes(col=model))+
      #geom_smooth(aes(fill=model), span=0.6)
      stat_summary(aes(col=model),fun=mean, geom="line")
    
    ggplot(data=df, aes(ang_S, ang_M, col=model))+
      geom_point(aes(col=model))+
      #geom_smooth(aes(fill=model), span=0.6)
      stat_summary(aes(col=model),fun=mean, geom="line")+
      geom_abline(color = "orange")
    
    
    ggplot(data=df, aes(siz_S, absdist_MS, col=model))+
      geom_point(aes(col=model))+
      #geom_smooth(aes(fill=model), span=0.6)
      stat_summary(aes(col=model),fun=mean, geom="line")

    ggplot(data=df, aes(siz_S, ang_M, col=model))+
      geom_point(aes(col=model))+
      #geom_smooth(aes(fill=model), span=0.6)
      stat_summary(aes(col=model),fun=mean, geom="line")+
      geom_point(aes(col=model))
      geom_hline(yintercept=-0.35, color = "orange")

      # ggplot(data=df, aes(ecc_S, absdist_MS, col=model))+
      #   geom_point(aes(col=model))+
      #   #geom_smooth(aes(fill=model), span=0.6)
      #   stat_summary(aes(col=model),fun=mean, geom="line")
      # 
      # ggplot(data=df, aes(ecc_S, ang_M, col=model))+
      #   geom_point(aes(col=model))+
      #   #geom_smooth(aes(fill=model), span=0.6)
      #   stat_summary(aes(col=model),fun=mean, geom="line")+
      #   geom_point(aes(col=model))+
      # geom_hline(yintercept=-0.35, color = "orange")
      
    # ggplot(data=df, aes(ecc_S, absdist_MS, col=model))+
    #   geom_point(aes(col=model))+
    #   geom_smooth(aes(fill=model), span=0.45)
    #   #stat_summary(aes(col=model),fun=mean, geom="line")    

dev.off()
