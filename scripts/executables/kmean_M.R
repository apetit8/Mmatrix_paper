source("../functions_R/network.R")
source("../functions_R/tools.R")
library(FactoMineR)
library(factoextra)
library(gridExtra)
#######################
sims.dirs <- list.dirs("../../simul/cercle/wag2/c1_1")
sims.dir  <- sims.dirs[2:30]
of        <- "sc1_1"
what      <- "angle" #"angle" or "ecc"
where     <- "cercle"
grp <- 6
scaling <- TRUE
#######################

df.topo <- df.topo.raw(sims.dir)


names(df.topo)[names(df.topo) == "V9"] <- "A-B"
names(df.topo)[names(df.topo) == "V10"] <- "A-C"
names(df.topo)[names(df.topo) == "V11"] <- "A-D"
names(df.topo)[names(df.topo) == "V12"] <- "B-A"
names(df.topo)[names(df.topo) == "V14"] <- "B-C"
names(df.topo)[names(df.topo) == "V15"] <- "B-D"
names(df.topo)[names(df.topo) == "V16"] <- "C-A"
names(df.topo)[names(df.topo) == "V17"] <- "C-B"
names(df.topo)[names(df.topo) == "V19"] <- "C-D"
names(df.topo)[names(df.topo) == "V20"] <- "D-A"
names(df.topo)[names(df.topo) == "V21"] <- "D-B"
names(df.topo)[names(df.topo) == "V22"] <- "D-C"

dfx <- df.topo[,c(9,10,11,12,14,15,16,17,19,20,21,22)]
dfy <- as.data.frame(df.topo$ang_M)
df <-as.data.frame(cbind(dfx,dfy))
names(df)[names(df) == "df.topo$ang_M"] <- "M_angle"
df[] <- lapply( df, as.numeric)
df.topo[] <- lapply( df.topo, as.numeric)

pdf(sprintf("../../figures/%s/%s_kmean_M.pdf", where, of), width=6, height=3)
#Determine number of clusters
wss <- (nrow(df[,1:12])-1)*sum(apply(df[,1:12],2,var))
for (i in 2:15) wss[i] <- sum(kmeans(df[,1:12],
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

fviz_nbclust(df[,1:12], kmeans, method = "wss") +
  geom_vline(xintercept = grp, linetype = 2)+
  labs(subtitle = "Elbow method")

# K-Means Cluster Analysis
fit <- kmeans(df[,1:12], grp) # grp cluster solution
# get cluster means
aggregate(df[,1:12],by=list(fit$cluster),FUN=mean)
# append cluster assignment
df <- data.frame(df, fit$cluster)
df$fit.cluster <- as.factor(df$fit.cluster)
df.topo <- data.frame(df.topo, df$fit.cluster)

ggplot(data=df.topo, aes(as.factor(round(x=ang_M, digits = 1)) , y=frequency(ang_M), Freq, fill = as.factor(df.fit.cluster)))+
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#02db77","#6ae8b1", "#8862ba", "#eb8698", "#E69F00", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90))+
  labs(title=sprintf("Clusters distribution in %s", of), x ="M angle", y = "Occurences /1160")
dev.off()


