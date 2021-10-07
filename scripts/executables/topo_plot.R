source("../functions_R/network.R")
source("../functions_R/tools.R")
library(ggplot2)
#######################
sims.dirs <- list.dirs("../../simul/Wagner/topos/exploration/angle_0.8_0.2")
sims.dir  <- sims.dirs[2:30]
of        <- "l4_0.8_0.2"
where     <- "topo"
treshold  <- 0.05
#######################


alltopos <- data.frame(Var1=character(), Feq=numeric())

for (i in sims.dir) {
  files     <- list.files(path = i, full.names=TRUE, pattern = "^simul")
  mytopos <- lapply(files, function(ff) {
    #print(ff)
    tt <- read.table(ff, header=TRUE)
    Wmat <- extract.W.matrix(tt)
    cWmat <- cleanW(Wmat, a=treshold)
    untopo <-unique.topo(cWmat)
  })
  topos <- data.frame(table(sapply(mytopos, function(topo) paste(topo, collapse="."))))
  alltopos <- rbind(alltopos, topos)
  print(sum(alltopos$Freq))
}

pdfname  <- print(sprintf("../../figures/%s/%s_topos_plot_%s.pdf", where, of, treshold))
pdf(pdfname, width=7, height=7)
layout(t(1:1))
  ggplot(as.data.frame(alltopos), aes(x=Var1, y=Freq, Freq, fill = Freq)) +     
    geom_col(position = 'dodge') + theme(axis.text.x = element_text(angle = 90)) + geom_bar(stat="identity", color="white") +
    labs(title=of, x ="Topology", y = "Effectif")
dev.off()

print("Alltopos done !")
