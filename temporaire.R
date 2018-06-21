### Here is the script standard to execute all your operations for one species from the plots extractions to the models running. 

rm(list = ls())
gc()
require(ade4)
library(sjPlot) 
library(lattice)
library(latticeExtra)
library(stringr)
library(data.table)
library(pscl)
library(MASS)
library(lme4)
library("glmmTMB")
library("bbmle")
library(piecewiseSEM)
library(pgirmess)
library(MuMIn)
library(spaMM)

# Extraction de la base de donnée climatique
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/6.climateSpecies.R"))
Allcode <- c("PICABI","PINPINA","FAGSYL","PINHAL","QUEROB","PINNIG","QUEPET","CASSAT","ABIALB","QUEPUB","QUEPYR","FRAEXC","PINPIN",
             "QUESUB","PINRAD","BETPEN","PINCAN","PINUNC","EUCGLO","QUEFAG","ALNGLU","POPTRE","ACEPSE","LARDEC","ROBPSE","POPNIG","PICSIT","QUERUB")
for (code in Allcode){
  try(ExtractClimate(CODE=code,yearINI=yearI,yearFIN=yearF,Intervalle=Inter,files.wanted=files.w,save=T),silent=T) # Maybe not all the intervalle
}

# Spatialisation of the wanted species 
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/8.spatialisation.R"))
Allcode <- c("FAGSYL","PINHAL","QUEROB","PINNIG","QUEPET","CASSAT","ABIALB","QUEPUB","QUEPYR","FRAEXC","PINPIN",
             "QUESUB","BETPEN","ALNGLU","POPTRE","ACEPSE","LARDEC","POPNIG")
Errors.files <- file(paste0(Dir,"species/Errors.spatialisation.Rout"), open="wt")
sink(file=Errors.files,append = T, type="message")
for (code in Allcode){
  try(Spatialisation(CODE = code,SaveAll = T,myPointsAll = T),silent=F) #re run the code with this arguments 
}
sink(type="message")
close(Errors.files)


# Save the Jpeg that did not work 
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/8.spatialisation.R"))
Allcode <- c("PINHAL","QUEROB","PINNIG","QUEPET","CASSAT","ABIALB","QUEPYR","FRAEXC","PINPIN",
             "QUESUB","BETPEN")
for (code in Allcode){
  i = 1
  for (i in 1:length(eumedclim.vars)) {
    assign(paste0(eumedclim.vars[i],"_T_Dist"),raster(list.files(paste0(Dir,"species/",code,""),pattern=paste0("^",eumedclim.vars[i],"_T_Dist.tif"), full.names=T)),envir = globalenv())
    try(Plot.map(CODE = code,Climvar = eumedclim.vars[i], Maptype=type[1], top = T, Save = T, myPoints = F),silent=F) #re run the code with this arguments 
  }
}
list.files(paste0(Dir,"species/",CODE,""),pattern=paste0(eumedclim.vars[1],"_T_Dist.tif"), full.names=T)


# Cut marginality for all species and do it again for the species whose plots number is too low. 
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/9.Marginality.R"))
Allcode <- list("PICABI","PINPINA","QUEILE","POPTRE","LARDEC","POPNIG")
Errors.files2 <- file(paste0(Dir,"species/Errors.Marginality.Rout"), open="wt")
sink(Errors.files2, type="message")
for (code in Allcode){
  try(Acp1000(CODE = code,nsample=10000,NF=2,seuil=0.8,seuilC=0.6,save=T),silent=T) #re run the code with this arguments # First function for PCA
}
Allcode <- c("PICABI","PINPINA","FAGSYL","PINHAL","QUEROB","QUEILE","PINNIG","QUEPET","CASSAT","ABIALB","QUEPUB","QUEPYR","FRAEXC","PINPIN",
             "QUESUB","BETPEN","ALNGLU","POPTRE","ACEPSE","LARDEC","POPNIG")
for (code in Allcode){
  try(Marginality_Levels(CODE=code,seuil=0.8,seuilC=0.6,save=T),silent=T)
}
sink(type="message")
close(Errors.files2)
mclapply(Allcode, function(x){Acp1000(CODE=x,nsample=10000,NF=2,seuil=0.8,seuilC=0.6,save=T)},mc.cores = 6)
mclapply(Allcode, function(x){Marginality_Levels(CODE = x,seuil=0.8,seuilC=0.6,save=T)},mc.cores = 4)

Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/10.Dfplots.R"))
Allcode <- list("LARDEC")
Errors.files <- file(paste0(Dir,"our-data/species/Errors.dfplots.Rout"), open="wt")
sink(Errors.files, type="message")
for (code in Allcode){
  try(MyDFs(dir="bureau",CODE = code,seuil=0.8,seuilC=0.6),silent=T) #re run the code with this arguments # First function for PCA
}
sink(type="message")
close(Errors.files)


Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/10bis.SPEI.R"))
Years <- c("spei01","spei03","spei06","spei18","spei24","spei36","spei48")
for (years in Years){
  try(SPEI(nStart = nStart,nEnd = nEnd, ncname = years),silent=T)
}



Dir = c("/Users/alexandrechangenet/Dropbox/FUNDIV/") # Directory 
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")

setwd(Dir)
# Load the function 
source(paste0(Dir,"Myscripts/function1.R.squared.r"))
source(paste0(Dir,"Myscripts/function2.Saving.R"))
source(paste0(Dir,"Myscripts/function3.ModelBoot.R"))
source(paste0(Dir,"Myscripts/function4.Diagnostic.R"))
source(paste0(Dir,"Myscripts/function5.Premodel.R"))

#Extraction of the database of the wanted species : 

CODE = "BETPEN" 
#"ALNGLU" 
#"ABIALB"
#"BETPEN"
#"PICABI"
#"PINPINA"
seuil = 0.80
#Dir = c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models")) # Directory 
Dir = paste0(Dir,"our-data/species/",CODE,"/CLIMAP") # Directory 
setwd(Dir)
list.files(Dir,pattern = paste0(seuil,".rds"))
dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
#dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
dir.create(path=paste0(Dir,"/Models"))
dir.create(path=paste0(Dir,"/Models/Negbin"))
dir.create(path=paste0(Dir,"/Models/binomial"))
Dir =c(paste0(Dir,"/Models"))
setwd(Dir)

Premodel(z=dfplot2,Resp=Resp,Explain=Explain,size=4,save=T)
Dir =c(paste0(Dir,"/binomial"))
setwd(Dir)
Mbin2BETPEN <- fitme(sp.mort.bin ~ treeNbr + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + bio14_climate_mean.30 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BAj.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:bio14_climate_mean.30 + BAj.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + (1|country), data=dfplot2,family = binomial(),method='REML')
Mbin80EtestOffsetPINSYL <- fitme(sp.mort.bin ~ BAIj.plot.bis + treeNbr +  offset(log(yearsbetweensurveys)) + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + bio14_climate_mean.30 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + (1|country), data=dfplot2, family = binomial(link = "cloglog"),method='REML')
Saving(Mbin2BETPEN)
Extraction(Mbin2BETPEN)
ModelBoot(Mbin80EtestOffsetPINSYL,5,7,LvL=50,CAT=CAT,nBoot=50,Yportion = 0.80,saveboot=T,nCoeur=10)
Diagnostic(Mbin80EtestOffsetPINSYL,0.66,F)

model <- Mbin1BETPEN
link <- model$family$link
family. <- model$family$family[1]

Dir =c(paste0(Dir,"/Negbin"))
setwd(Dir)
MnbZT1BETPEN <- fitme(sp.mortality.plot.count.yr ~ treeNbr + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + tmean.djf_climate_mean.30 + ppet.mean_climate_mean.30 + Plotcat + I(tmean.djf_climate_mean.30^2) + I(ppet.mean_climate_mean.30^2) + BA.ha.plot.1:tmean.djf_climate_mean.30 + BAj.plot.1:tmean.djf_climate_mean.30 + BA.ha.plot.1:ppet.mean_climate_mean.30 + BAj.plot.1:ppet.mean_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + tmean.djf_climate_mean.30:ppet.mean_climate_mean.30 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:tmean.djf_climate_mean.30 + Plotcat:ppet.mean_climate_mean.30 + (1|country), data=subset(dfplot2,sp.mortality.plot.count.yr>0),family=negbin(),method='REML')
MnbZT.Ibis.PINSYL <- fitme(sp.mortality.plot.count.yr ~ BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + tmean.djf_climate_mean.30 + ppet.mean_climate_mean.30 + Plotcat + I(ppet.mean_climate_mean.30^2) + BAj.plot.1:tmean.djf_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + (1|country), data=subset(dfplot2,sp.mortality.plot.count.yr>0),family=negbin(),method='REML')
Saving(MnbZT.Ibis.PINSYL)
Extraction(MnbZT.Ibis.PINSYL)
ModelBoot(MnbZT.Ibis.PINSYL,3,6,LvL=20,CAT,nBoot=30,Yportion = 0.80,saveboot=F,nCoeur=10)
Diagnostic(MnbZT.Ibis.PINSYL,0.66,F)




list.dirs(getwd())
setwd(list.dirs(getwd())[27]) # Chose the model file
Mbin80EtestOffset <- get(load(file = "Mbin80EtestOffset.rda")) 
rm("x")
Diagnostic(Mbin80EtestOffset,0.66,T)
# bin
dfplot2=rbind(sample_n(dfplot2[dfplot2$sp.mort.bin==0 & dfplot2$Plotcat.80.60==0,],500),
              sample_n(dfplot2[dfplot2$sp.mort.bin==0 & dfplot2$Plotcat.80.60==1,],500),
              sample_n(dfplot2[dfplot2$sp.mort.bin==0 & dfplot2$Plotcat.80.60==2,],150),
              sample_n(dfplot2[dfplot2$sp.mort.bin==1 & dfplot2$Plotcat.80.60==0,],500),
              sample_n(dfplot2[dfplot2$sp.mort.bin==1 & dfplot2$Plotcat.80.60==1,],500),
              sample_n(dfplot2[dfplot2$sp.mort.bin==1 & dfplot2$Plotcat.80.60==2,],51))

# Negbin 
dfplot2=rbind(sample_n(dfplot2[dfplot2$sp.mortality.plot.count.yr>0 & !is.na(dfplot2$sp.mortality.plot.count.yr) & dfplot2$Plotcat.80.60==0,],400),
              sample_n(dfplot2[dfplot2$sp.mortality.plot.count.yr>0 & !is.na(dfplot2$sp.mortality.plot.count.yr) & dfplot2$Plotcat.80.60==1,],400),
              sample_n(dfplot2[dfplot2$sp.mortality.plot.count.yr>0 & !is.na(dfplot2$sp.mortality.plot.count.yr) & dfplot2$Plotcat.80.60==2,],40))


dfplot2$sp.mortality.plot.ba <- round((dfplot2$sp.mortality.plot.ba*1000/dfplot2$yearsbetweensurveys))
MnbZT80.60.M.BA <- fitme(sp.mortality.plot.ba ~ BAIj.plot.bis +  treeNbr + BA.ha.plot.1 + BAj.plot.1 + tmean.djf_climate_mean.30 + Plotcat.80.60 + I(tmean.djf_climate_mean.30^2) + I(ppet.mean_climate_mean.30^2) + BAj.plot.1:tmean.djf_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BA.ha.plot.1:Plotcat.80.60 + BAj.plot.1:Plotcat.80.60 + Plotcat.80.60:tmean.djf_climate_mean.30 + (1|country), data=subset(dfplot2,sp.mortality.plot.ba>0),family=negbin(),method='REML')
Saving(MnbZT80.60.M.BA)
MnbZT80.60.M.BAcontrol <- fitme(sp.mortality.plot.ba ~ BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + tmean.djf_climate_mean.30 + Plotcat.80.60 + I(tmean.djf_climate_mean.30^2) + I(ppet.mean_climate_mean.30^2) + BAj.plot.1:tmean.djf_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BA.ha.plot.1:Plotcat.80.60 + BAj.plot.1:Plotcat.80.60 + Plotcat.80.60:tmean.djf_climate_mean.30 + (1|country), data=subset(dfplot2,sp.mortality.plot.ba>0),family=negbin(),method='REML')
Saving(MnbZT80.60.M.BAcontrol)
MnbZT80.60.M <- fitme(sp.mortality.plot.count ~ BAIj.plot.bis +  offset(log(yearsbetweensurveys)) + treeNbr + BA.ha.plot.1 + BAj.plot.1 + tmean.djf_climate_mean.30 + Plotcat.80.60 + I(tmean.djf_climate_mean.30^2) + I(ppet.mean_climate_mean.30^2) + BAj.plot.1:tmean.djf_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BA.ha.plot.1:Plotcat.80.60 + BAj.plot.1:Plotcat.80.60 + Plotcat.80.60:tmean.djf_climate_mean.30 + (1|country), data=subset(dfplot2,sp.mortality.plot.rate.yr>0),family=negbin(),method='REML')
Saving(MnbZT80.60.M)
MnbZT80.60.Mbis <- fitme(sp.mortality.plot.rate.yr ~ BAIj.plot.bis + treeNbr + BA.ha.plot.1 + BAj.plot.1 + tmean.djf_climate_mean.30 + Plotcat.80.60 + I(tmean.djf_climate_mean.30^2) + I(ppet.mean_climate_mean.30^2) + BAj.plot.1:tmean.djf_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BA.ha.plot.1:Plotcat.80.60 + BAj.plot.1:Plotcat.80.60 + Plotcat.80.60:tmean.djf_climate_mean.30 + (1|country), data=subset(dfplot2,sp.mortality.plot.rate.yr>0),family=negbin(),method='REML')
Saving(MnbZT80.60.Mbis)
MnbZT80.60.Mbiscontrol <- fitme(sp.mortality.plot.rate.yr ~ BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + tmean.djf_climate_mean.30 + Plotcat.80.60 + I(tmean.djf_climate_mean.30^2) + I(ppet.mean_climate_mean.30^2) + BAj.plot.1:tmean.djf_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BA.ha.plot.1:Plotcat.80.60 + BAj.plot.1:Plotcat.80.60 + Plotcat.80.60:tmean.djf_climate_mean.30 + (1|country), data=subset(dfplot2,sp.mortality.plot.rate.yr>0),family=negbin(),method='REML')
Saving(MnbZT80.60.Mbiscontrol)

spaMM.options(separation_max=2200)
Mbin80EtestOffset <- fitme(sp.mort.bin ~ BAIj.plot.bis + treeNbr +  offset(log(yearsbetweensurveys)) + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + bio14_climate_mean.30 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + (1|country), data=dfplot2, family = binomial(link = "cloglog"),method='REML')
Saving(Mbin80EtestOffset)
Mbin80Econtrol <- fitme(sp.mort.bin ~ BAIj.plot.bis + treeNbr + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + bio14_climate_mean.30 + Plotcat.80.60 + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BAj.plot.1:Plotcat.80.60 + Plotcat.80.60:bio1_climate_mean.30 + Plotcat.80.60:bio14_climate_mean.30 + (1|country), data=dfplot2, family = binomial,method='REML')
Saving(Mbin80Econtrol)
Mbin80Econtrol2 <- fitme(sp.mort.bin ~ BAIj.plot.bis + treeNbr + yearsbetweensurveys + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + bio14_climate_mean.30 + Plotcat.80.60 + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BAj.plot.1:Plotcat.80.60 + Plotcat.80.60:bio1_climate_mean.30 + Plotcat.80.60:bio14_climate_mean.30 + (1|country), data=dfplot2, family = binomial,method='REML')
Saving(Mbin80Econtrol2)


list.dirs(getwd())
setwd(list.dirs(getwd())[27]) # Chose the model file
Mbin80EtestOffset <- get(load(file = "Mbin80EtestOffset.rda")) 
rm("x")
y <- deparse(substitute(Mbin80EtestOffset))
Extraction(get(y))
ModelBoot(Mbin80EtestOffset,5,8,20,CAT=CAT,30,T)
Diagnostic(Mbin80EtestOffset,0.66,T)


