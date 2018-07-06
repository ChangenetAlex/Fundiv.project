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


###################################################################################
###                                                                            ####
###                                                                            #### # Here, 6 scripts are run
###     PART two : how to extract the information we want for all species      ####
###                                                                            ####
###                                                                            ####
###################################################################################
Dir = c("/Users/alexandrechangenet/Dropbox/FUNDIV/") # Directory 
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")

setwd(Dir)
# Load the function 
source(paste0(Dir,"Myscripts/Fundiv.project/function1.R.squared.r"))
source(paste0(Dir,"Myscripts/Fundiv.project/function2.Saving.R"))
source(paste0(Dir,"Myscripts/Fundiv.project/function3.ModelBoot.R"))
source(paste0(Dir,"Myscripts/Fundiv.project/function4.Diagnostic.R"))
source(paste0(Dir,"Myscripts/Fundiv.project/function5.Premodel.R"))
source(paste0(Dir,"Myscripts/Fundiv.project/function6.Coefs.R"))

#Extraction of the database of the wanted species : 

CODE = "FAGSYL" 
#"ALNGLU" 
#"ABIALB"
#"BETPEN"
#"PICABI"
#"PINPINA"
#"FAGSYL"
#"PINHAL"
#"QUEROB"
#"QUEILE"
#"PINNIG"
#"QUEPET"
#"CASSAT"
#"ABIALB"
#"QUEPUB"
#"QUEPYR"
#"FRAEXC"
#"PINPIN"
#"QUESUB"
#"BETPEN"
#"ALNGLU"
#"POPTRE"
#"ACEPSE"
#"LARDEC"
#"POPNIG"

# Load the df

seuil = 0.70
#Dir = c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models")) # Directory 
Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
setwd(Dir)
list.files(Dir,pattern = paste0(seuil,".rds"))
dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
#dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
#dir.create(path=paste0(Dir,"/Models"))
#dir.create(path=paste0(Dir,"/Models/Negbin"))
#dir.create(path=paste0(Dir,"/Models/binomial"))
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models"))
setwd(Dir)



# Perform the anaysis on all variable three by three 
for (i in 1:((length(Explain)-1)/3)){
  Explain <- c("BA.ha.plot.1","BA.O.plot.1","BAj.plot.1","dbh.plot.mean","BAIj.plot.1.mean","BAIj.plot.1",
               "logBA.ha.plot.1","logBA.O.plot.1","logBAj.plot.1","logdbh.plot.mean","logBAIj.plot.1.mean","logBAIj.plot.1",
               "sqrtBA.ha.plot.1","sqrtBA.O.plot.1","sqrtBAj.plot.1","sqrtdbh.plot.mean","sqrtBAIj.plot.1.mean","sqrtBAIj.plot.1",
               "treeNbr","yearsbetweensurveys","bio1_climate_mean.30","bio14_climate_mean.30","min_spei12","mean_spei12","Plotcat")
  Explain <- Explain[c((3*i-2):(3*i),length(Explain))]
  Premodel(z=dfplotbis,Resp=Resp,Explain=Explain,size=5,save=T) #
  Premodel(z=dfplot2,Resp=Resp,Explain=Explain,size=5,save=T)   # Apply the functiçon on both database on all trasnformed and no trasnformed data
} 
# Premodel analysis for 12 variables and to obtain the global VIF indices
Explain <- c("bio1_climate_mean.30","bio14_climate_mean.30","min_spei12","mean_spei12",
             "BA.ha.plot.1","BA.O.plot.1","BAj.plot.1","dbh.plot.mean","BAIj.plot.1.mean","BAIj.plot.1","treeNbr","yearsbetweensurveys","Plotcat")
Resp <- c("sp.mort.bin")
#Resp <- c("sp.mortality.plot.count.yr")
Premodel(z=dfplot2,Resp=Resp,Explain=Explain,size=3,save=T) # Obtain the global VIF for all variable (not transformed)
  



Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/binomial"))
setwd(Dir)


hist(dfplotbis$BA.ha.plot.1,breaks=100)
hist(dfplotbis$logBA.ha.plot.1,breaks=100)
hist(dfplotbis$sqrtBA.ha.plot.1,breaks=100)
hist(dfplot2$BA.ha.plot.1,breaks=100)
hist(dfplot2$logBA.ha.plot.1,breaks=100)
hist(dfplot2$sqrtBA.ha.plot.1,breaks=100) # To test 

hist(dfplotbis$BA.O.plot.1,breaks=100)
hist(dfplotbis$logBA.O.plot.1,breaks=100)
hist(dfplotbis$sqrtBA.O.plot.1,breaks=100)
hist(dfplot2$BA.O.plot.1,breaks=30)
hist(dfplot2$logBA.O.plot.1,breaks=30)
hist(dfplot2$sqrtBA.O.plot.1,breaks=30) # Mieux mais zero inflated 

hist(dfplotbis$BAj.plot.1,breaks=50)
hist(dfplotbis$logBAj.plot.1,breaks=50)
hist(dfplotbis$sqrtBAj.plot.1,breaks=50)
hist(dfplot2$BAj.plot.1,breaks=50)
hist(dfplot2$logBAj.plot.1,breaks=100) # Mieux mais zero inflated
hist(dfplot2$sqrtBAj.plot.1,breaks=100)


hist(dfplotbis$BAIj.plot.1.mean,breaks=50)
hist(dfplotbis$logBAIj.plot.1.mean,breaks=50)
hist(dfplotbis$sqrtBAIj.plot.1.mean,breaks=50)
hist(dfplot2$BAIj.plot.1.mean,breaks=50)
hist(dfplot2$logBAIj.plot.1.mean,breaks=50)
hist(dfplot2$sqrtBAIj.plot.1.mean,breaks=50)  ## Best 


hist(dfplotbis$BAIj.plot.1,breaks=50)
hist(dfplotbis$logBAIj.plot.1,breaks=50)
hist(dfplotbis$sqrtBAIj.plot.1,breaks=50)
hist(dfplot2$BAIj.plot.1,breaks=50)
hist(dfplot2$logBAIj.plot.1,breaks=50)
hist(dfplot2$sqrtBAIj.plot.1,breaks=50) # Best 


hist(dfplotbis$dbh.plot.mean,breaks=50)
hist(dfplotbis$logdbh.plot.mean,breaks=50)
hist(dfplotbis$sqrtdbh.plot.mean,breaks=50)
hist(dfplot2$dbh.plot.mean,breaks=50)
hist(dfplot2$logdbh.plot.mean,breaks=30) # Best 
hist(dfplot2$sqrtdbh.plot.mean,breaks=30)

# Models 
Mbin1FAGSYL <- fitme(sp.mort.bin ~ treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1 + bio14_climate_mean.30 + logbio14 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1 + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:logbio14 + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:logbio14 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1:logbio14 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:logbio1 + Plotcat:logbio14 + (1|country), data=dfplot2,family = binomial(),method='REML')
Mbin1FAGSYL <- fitme(sp.mort.bin ~ treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1 + bio14_climate_mean.30 + logbio14 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1 + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:logbio14 + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:logbio14 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1:logbio14 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:logbio1 + Plotcat:logbio14 + (1|country), data=dfplot2,family = binomial(),method='REML')

Mbin2FAGSYLspaMM <- fitme(sp.mort.bin ~ treeNbr + yearsbetweensurveys + min_spei12 + mean_spei12 + min_spei24 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1_climate_mean.30 + bio14_climate_mean.30 + logbio14_climate_mean.30 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1_climate_mean.30 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1_climate_mean.30 + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:logbio14_climate_mean.30 + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:logbio14_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1_climate_mean.30:logbio14_climate_mean.30 + BA.ha.plot.1:min_spei12 + BA.ha.plot.1:mean_spei12 + BAj.plot.1:min_spei12 + BAj.plot.1:mean_spei12 + bio1_climate_mean.30:min_spei12 + bio1_climate_mean.30:mean_spei12 + bio14_climate_mean.30:min_spei12 + bio14_climate_mean.30:mean_spei12 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:logbio1_climate_mean.30 + Plotcat:logbio14_climate_mean.30 + Plotcat:min_spei12 + Plotcat:mean_spei12 + (1 | country), data=dfplot2,family = binomial(),method='REML')
Mbin3FAGSYLspaMM <- fitme(sp.mort.bin ~ treeNbr + yearsbetweensurveys + min_spei12 + mean_spei12 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1_climate_mean.30 + bio14_climate_mean.30 + logbio14_climate_mean.30 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1_climate_mean.30 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1_climate_mean.30 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1_climate_mean.30:logbio14_climate_mean.30 + BA.ha.plot.1:min_spei12 + BA.ha.plot.1:mean_spei12 + bio1_climate_mean.30:min_spei12 + BA.ha.plot.1:Plotcat + Plotcat:min_spei12 + Plotcat:mean_spei12 + (1 | country), data=dfplot2,family = binomial(),method='REML')

Mbin1FAGSYL <- fitme(sp.mort.bin ~ treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + log(bio1_climate_mean.30+5) + bio14_climate_mean.30 + log(bio14_climate_mean.30+5) + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:log(bio1_climate_mean.30+5) + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:log(bio1_climate_mean.30+5) + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:log(bio14_climate_mean.30+5) + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:log(bio14_climate_mean.30+5) + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + log(bio1_climate_mean.30+5):log(bio14_climate_mean.30+5) + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:log(bio1_climate_mean.30+5) + Plotcat:log(bio14_climate_mean.30+5) + (1|country), data=dfplot2,family = binomial(),method='REML')
Mbin1FAGSYL <- fitme(sp.mort.bin ~ treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + bio14_climate_mean.30 + log(dfplot2$bio14_climate_mean.30+5) + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:log(bio1_climate_mean.30+1) + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:log(bio1_climate_mean.30+1) + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:log(bio14_climate_mean.30+1) + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:log(bio14_climate_mean.30+1) + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + log(bio1_climate_mean.30+1):log(bio14_climate_mean.30+1) + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:log(bio1_climate_mean.30+1) + Plotcat:log(bio14_climate_mean.30+1) + (1|country), data=dfplot2,family = binomial(),method='REML')
Mbin80EtestOffsetPINSYL <- fitme(sp.mort.bin ~ BAIj.plot.bis + treeNbr +  offset(log(yearsbetweensurveys)) + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + bio14_climate_mean.30 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + (1|country), data=dfplot2, family = binomial(link = "cloglog"),method='REML')

### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[2]) # Chose the model file
list.files()
Mbin31 <- get(load(file = "Mbin31.rda")) 
rm("x")


# Process the models : bootstrap and coefficient, saving etc 
Saving(Mbin3FAGSYLspaMM) # Evaluation of the model 
Extraction(Mbin31)
ModelBoot(Mbin3FAGSYLspaMM,4,12,LvL=30,CAT=CAT,nBoot=10,Yportion = 0.80,saveboot=T,nCoeur=10)
Diagnostic(Mbin3FAGSYLspaMM,0.66,F)
ExtracTest(Mbin31) # Me donne la liste des paramètres de mon modèles. 
for (i in c("BAIj.plot.1","BAIj.plot.1.mean","dbh.plot.mean","BA.ha.plot.1","BAj.plot.1","treeNbr","bio1_climate_mean.30","bio12_climate_mean.30","min_spei12","mean_spei12")){
  Effect_coef(Mbin31,i)}     # Para que l'on veut regarder  parmis ceux citer  
Effect_summary(Mbin31,"competition") #
ggEffect(Mbin31,'ABS',"indiv",band=T)
