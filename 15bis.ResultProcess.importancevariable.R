library(reshape2)
rm(list = ls())
gc()
getwd()

tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresults.csv")
summary(as.factor(tableresults$variable))
#tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
tableresults <- tableresults[-c(which(tableresults$variable=="(Intercept)")),] # On enlève l'intercept
#tableresults$signif <- as.character(tableresults$signif)
tableresults$species <- paste0(tableresults$eco," ",tableresults$species)
tableresults <- tableresults[,c(1,2,3)] # remove conif and angio
tableresults$estimates <- round(tableresults$estimates,3)
#tableresults$eco <- as.character(tableresults$eco)
tableresults$variable <- as.character(tableresults$variable)
tableresults$variable <- sub(":Plotcat1"," : Marginality (TE)",tableresults$variable)
tableresults$variable <- sub("Plotcat1:","Marginality (TE) : ",tableresults$variable)
tableresults$variable <- sub(":Plotcat2"," : Marginality (LE)",tableresults$variable)
tableresults$variable <- sub("Plotcat2:","Marginality (LE) : ",tableresults$variable)
tableresults$variable <- sub("Plotcat2","Marginality (LE)",tableresults$variable)
tableresults$variable <- sub("Plotcat1","Marginality (TE)",tableresults$variable)
tableresults$variable <- sub("_climate_mean.30","",tableresults$variable)
tableresults$variable <- sub("bio12","Precipitation1",tableresults$variable)
tableresults$variable <- sub("bio13","Precipitation2",tableresults$variable)
tableresults$variable <- sub("bio14","Precipitation3",tableresults$variable)
tableresults$variable <- sub("bio1","Temperature1",tableresults$variable)
tableresults$variable <- sub("bio5","Temperature2",tableresults$variable)
tableresults$variable <- sub("mean_spei12","Drought1",tableresults$variable)
tableresults$variable <- sub("min_spei12","Drought2",tableresults$variable)
tableresults$variable <- sub("ppet.mean","Drought3",tableresults$variable)
tableresults$variable <- sub("sqrtBA.ha.plot.1","Competition1",tableresults$variable)
tableresults$variable <- sub("sqrtBA.O.plot.1","Competition2",tableresults$variable)
tableresults$variable <- sub("logBAj.plot.1","Competition3",tableresults$variable)
tableresults$variable <- sub("log","",tableresults$variable)
tableresults$variable <- sub("sqrt","",tableresults$variable)

datDF=split(tableresults[,],as.character(tableresults$species)) # split mon df
A=lapply(datDF, function(x) x[order(abs(x[,3]),decreasing=T),][1:5,])
A <- do.call(rbind,A)
A$estimates <- paste0(A$estimates," ",A$variable)
A <- A[,c(1,3)]
A$Varank <- c("Estimate 1","Estimate 2","Estimate 3","Estimate 4","Estimate 5")
A <- dcast(A,species~Varank,value.var="estimates")
A
write.table(A, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/Coef.binom.csv",row.names = F)


#### Idem ZTNB models 
tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresultsNB.csv")
summary(as.factor(tableresults$variable))
#tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
tableresults <- tableresults[-c(which(tableresults$variable=="(Intercept)")),] # On enlève l'intercept
#tableresults$signif <- as.character(tableresults$signif)
tableresults$species <- paste0(tableresults$eco," ",tableresults$species)
tableresults <- tableresults[,c(1,2,3)] # remove conif and angio
tableresults$estimates <- round(tableresults$estimates,3)
#tableresults$eco <- as.character(tableresults$eco)
tableresults$variable <- as.character(tableresults$variable)
tableresults$variable <- sub(":Plotcat1"," : Marginality (TE)",tableresults$variable)
tableresults$variable <- sub("Plotcat1:","Marginality (TE) : ",tableresults$variable)
tableresults$variable <- sub(":Plotcat2"," : Marginality (LE)",tableresults$variable)
tableresults$variable <- sub("Plotcat2:","Marginality (LE) : ",tableresults$variable)
tableresults$variable <- sub("Plotcat2","Marginality (LE)",tableresults$variable)
tableresults$variable <- sub("Plotcat1","Marginality (TE)",tableresults$variable)
tableresults$variable <- sub("_climate_mean.30","",tableresults$variable)
tableresults$variable <- sub("bio12","Precipitation1",tableresults$variable)
tableresults$variable <- sub("bio13","Precipitation2",tableresults$variable)
tableresults$variable <- sub("bio14","Precipitation3",tableresults$variable)
tableresults$variable <- sub("bio1","Temperature1",tableresults$variable)
tableresults$variable <- sub("bio5","Temperature2",tableresults$variable)
tableresults$variable <- sub("tmean.djf","Temperature3",tableresults$variable)
tableresults$variable <- sub("mean_spei12","Drought1",tableresults$variable)
tableresults$variable <- sub("min_spei12","Drought2",tableresults$variable)
tableresults$variable <- sub("ppet.mean","Drought3",tableresults$variable)
tableresults$variable <- sub("sqrtBA.ha.plot.1","Competition1",tableresults$variable)
tableresults$variable <- sub("sqrtBA.O.plot.1","Competition2",tableresults$variable)
tableresults$variable <- sub("logBAj.plot.1","Competition3",tableresults$variable)
tableresults$variable <- sub("log","",tableresults$variable)
tableresults$variable <- sub("sqrt","",tableresults$variable)

datDF=split(tableresults[,],as.character(tableresults$species)) # split mon df
B=lapply(datDF, function(x) x[order(abs(x[,3]),decreasing=T),][1:5,])
B <- do.call(rbind,B)
B$estimates <- paste0(B$estimates," ",B$variable)
B <- B[,c(1,3)]
B$Varank <- c("Estimate 1","Estimate 2","Estimate 3","Estimate 4","Estimate 5")
B <- dcast(B,species~Varank,value.var="estimates")
B

write.table(B, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/Coef.negbin.csv",row.names = F)



tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/1_table_bin_gene.csv")
tableresults <- tableresults[,-c(6)]
tableresults[3:7] <- round(tableresults[,3:7],2)
write.table(tableresults, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/binomallmodels.csv",row.names = F)


tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/3_table_negbin_gene.csv")
tableresults <- tableresults[,-c(7)]
tableresults[,3:8] <- round(tableresults[,3:8],2)
write.table(tableresults, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/negbinmodels.csv",row.names = F)






