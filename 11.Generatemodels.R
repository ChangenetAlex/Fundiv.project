#########################################
#####   Work and Process data     #######
#########################################


metric = c("mean.","min.","max.") #mean
Interv = c(30) #5 10 15 30

#bioclim=paste("bio",c(1,12,13,14,2,5,6),"_climate_",sep="") #
#T.seas=paste0("tmean.",c("djf","jja","mam","son"),"_climate_")#
#P.seas=paste0("prec.",c("djf","jja","mam","son"),"_climate_")#
#pet=paste0("pet.",c("max","mean","min"),"_climate_") #
#ppet=paste0("ppet.",c("max","mean","min"),"_climate_")#
#eumedclim.vars=(c(bioclim, pet, ppet, P.seas, T.seas))
#Myvariable <- paste0(rep(eumedclim.vars,length(Interv)*length(metric)),rep(metric,length(Interv)*length(eumedclim.vars)),rep(Interv,length(eumedclim.vars)*length(metric)))

Allvariable.axe1 = c(paste0("bio1_climate_mean.",Interv), #Axe 1
                     paste0("tmean.djf_climate_mean.",Interv),
                     paste0("bio5_climate_mean.",Interv),
                     paste0("tmean.jja_climate_mean.",Interv))
Allvariable.axe2 = c(paste0("bio12_climate_mean.",Interv), #Axe 2
                     paste0("bio13_climate_mean.",Interv), #Axe 2
                     paste0("ppet.mean_climate_mean.",Interv),#Axe 2
                     paste0("bio14_climate_mean.",Interv))#Axe2

#Fonction pour écrire mes modèles : 
expand.grid(Allvariable.axe1,Allvariable.axe2)
n <- expand.grid(Allvariable.axe1,Allvariable.axe2,stringsAsFactors = F)
n <- n[n$Var1!=n$Var2,]

############################################
####                                    ####
#### Fonction to generate all the model ####
####       I want to evaluate :         ####
####                                    ####
############################################
i=1
vardep = c("sp.mortality.plot.rate.yr","sp.mort.bin") #keep this one 
varcat = c("Plotcat")
random <- c("country") #This will be kept and modified
for (i in 1:nrow(n)){
 varAll = c(n[i,1],n[i,2],"min_spei12","mean_spei12","BA.ha.plot.1","BAj.plot.1","BA.O.plot.1")
 vardep1 <- paste0("Mbin",i," <- fitme(",vardep[2]," ~ ",collapse = "")
 vardep2 <- paste0("MnbZT",i," <- fitme(",vardep[1]," ~ ",collapse = "")
 varAll1 <- paste0(varAll,collapse=" + ")
 varcat1 <- paste0(" + ",varcat,collapse="")
 carre <- paste0(" + I(",varAll[1:4],"^2)",collapse = "")
 interac1 <- expand.grid(varAll,varAll,stringsAsFactors = F)
 interac1 <- interac1[interac1$Var1!=interac1$Var2,]
 interac1 = paste0(" + ",interac1[,1],":",interac1[,2],collapse="")
 interac2 = paste0(" + ",expand.grid(varAll,varcat)[,1],":",expand.grid(varAll,varcat)[,2],collapse="")
 random1 <- paste0(" + (1|",random,")",collapse = "")
 finfun <- (", data=dfplot2, family = binomial,method='REML')")
 finfun1 <- (", data=subset(dfplot2,sp.mortality.plot.rate.yr>0),family=negbin(),method='REML')")
 capture.output(cat(paste0("#Mbin",i," = ",varAll[1]," ~ ",varAll[2]),"\n",paste0(vardep1,"BAIj.plot.bis + treeNbr + yearsbetweensurveys + ",varAll1," + ",varcat,carre,interac1," + ",interac2,random1,finfun,collapse = ""),"\n",sep=""),file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
 capture.output(cat(paste0("try(Saving(Mbin",i,"),silent=T)"),"\n","\n"),file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
 capture.output(cat(paste0("#MnbZT",i,"nb = ",varAll[1]," ~ ",varAll[2]),"\n",paste0(vardep2,"BAIj.plot.bis + treeNbr + yearsbetweensurveys + ",varAll1," + ",varcat,carre,interac1," + ",interac2,random1,finfun1,collapse = ""),"\n",sep=""),file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
 capture.output(cat(paste0("try(Saving(MnbZT",i,"),silent=T)"),"\n","\n"),file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
 i = i+1
}

