source("../functions_R/network.R")
source("../functions_R/tools.R")
library(ggplot2)
library(plyr)
library(tidyr)
#######################
sims.dirs <- list.dirs("../../simul/Wagner/scale0-1/ecc0.12/l4_0.5_0.8")
sims.dir  <- sims.dirs[2:8]
of        <- "l4_0.5_0.8"
what      <- "angle" #"angle" or "ecc"
where     <- "topo/ecc1"
threshold  <- 0.05
c <- 0.5 #constitutive expression
#######################

df.topo <- df.feat.and.topo(sims.dir, of, what, where, threshold, c)
#df.topo <- dftopo0.5

#palette
nb.cols <- length(unique(df.topo$topo))
mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)

pdfname  <- print(sprintf("../../figures/%s/%s_topos_feat_%s.pdf", where, of, threshold))
pdf(pdfname, width=20, height=10)
layout((1:2))
      ggplot(as.data.frame(df.topo), aes(x=topo, y=frequency(topo), Freq, fill = as.factor(df.topo$Value))) +
        theme(axis.text.x = element_text(angle = 90)) +
    geom_bar(stat="identity") +
        labs(title=of, x ="Topology", y = "Effectif")
    
    
    ggplot(as.data.frame(df.topo), aes(x=as.factor(df.topo$Value), y=frequency(df.topo$Value),  Freq, fill = topo)) +
      scale_fill_manual(values = mycolors)+
      theme(axis.text.x = element_text(angle = 90)) +
      geom_bar(stat="identity") +
      labs(title=of, x =what, y = "Effectif")
dev.off()



 #Occurrences table
 csvname   <- print(sprintf("../../figures/%s/%s_topos_list_%s.csv", where, of, threshold))
 occurences <- data.frame(sort(table(unlist(df.topo$topo)), decreasing = TRUE))
 write.table(data.frame(occurences), file = csvname, col.names = TRUE)
 sum(occurences$Freq)