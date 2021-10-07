###########################################
#import data and create a data frame containing only the 5000th generations
temp = list.files(pattern="*.txt")
simuls = lapply(temp, function(x)read.table(x, header=T)) #import list
simul = do.call("rbind", simuls) #assuming the same header/columns for all files
simul5000 = dplyr::filter(simul, Gen == 5000)
############################################


############################################
#Liste des correlations
corrG= (simul5000$CovPhen1_2)/sqrt((simul5000$VPhen2)*(simul5000$VPhen1))
plot(corrG)
boxplot(corrG) #Boxplot correlations

corrGmean = mean(corrG) #Correlations mean
