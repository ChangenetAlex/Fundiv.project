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

CODE = "ALNGLU" 
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

seuil = 0.60
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

###########################
## Premodel analysis ###### 
###########################

# Perform the anaysis on all variable three by three
#Resp <- c("sp.mortality.plot.count.yr")
Resp <- c("sp.mort.bin")
Explain <- c("BA.ha.plot.1","BA.O.plot.1","BAj.plot.1","dbh.plot.mean","BAIj.plot.1.mean","BAIj.plot.1",
             "logBA.ha.plot.1","logBA.O.plot.1","logBAj.plot.1","logdbh.plot.mean","logBAIj.plot.1.mean","logBAIj.plot.1",
             "sqrtBA.ha.plot.1","sqrtBA.O.plot.1","sqrtBAj.plot.1","sqrtdbh.plot.mean","sqrtBAIj.plot.1.mean","sqrtBAIj.plot.1",
             "treeNbr","yearsbetweensurveys","bio14_climate_mean.30","bio1_climate_mean.30","min_spei12","mean_spei12","Plotcat")

for (i in 1:((length(Explain)-1)/3)){
  Explain <- c("BA.ha.plot.1","BA.O.plot.1","BAj.plot.1","dbh.plot.mean","BAIj.plot.1.mean","BAIj.plot.1",
               "logBA.ha.plot.1","logBA.O.plot.1","logBAj.plot.1","logdbh.plot.mean","logBAIj.plot.1.mean","logBAIj.plot.1",
               "sqrtBA.ha.plot.1","sqrtBA.O.plot.1","sqrtBAj.plot.1","sqrtdbh.plot.mean","sqrtBAIj.plot.1.mean","sqrtBAIj.plot.1",
               "treeNbr","yearsbetweensurveys","bio14_climate_mean.30","bio1_climate_mean.30","min_spei12","mean_spei12","Plotcat")
  Explain <- Explain[c((3*i-2):(3*i),length(Explain))]
  print(Explain)
  Premodel(z=dfplot2,Resp=Resp,Explain=Explain,size=2,save=T)   # Apply the functictn on both database on all trasnformed and no trasnformed data
} 

# Premodel analysis for 12 variables and to obtain the global VIF indices
Explain <- c("bio1_climate_mean.30","bio14_climate_mean.30","min_spei12","mean_spei12",
             "sqrtBA.ha.plot.1","sqrtBA.O.plot.1","logdbh.plot.mean","sqrtBAIj.plot.1.mean","sqrtBAIj.plot.1","treeNbr","yearsbetweensurveys","Plotcat")
Premodel(z=dfplot2,Resp=Resp,Explain=Explain,size=3,save=T) # Obtain the global VIF for all variable (not transformed)


#########################
#####  Models Full  #####
#########################
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/binomial"))
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/Negbin"))

setwd(Dir)
### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[59]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[2], ignore.case = FALSE,fixed = T),get(load(file = list.files()[2]))) 
rm("x")
summary(get(sub(".rda","",list.files()[2])))


##########################
#### Coef + valid ########
##########################
Diagnostic(Mbin13A,0.66,F) # Pbm with the coreelog function 
ExtracTest(Mbin13A) # Me donne la liste des paramètres de mon modèles. 
for (i in c(A[c(1:11)])){
  Effect_coef(Mbin13A,i)}     # Para que l'on veut regarder  parmis ceux citer  
Effect_summary(Mbin13A,"biotic") #
ggEffect(Mbin13A,'REL',"sum",band=T)
ggEffect(Mbin13A,'ABS',"sum",band=T) 

###########################
#####    Affined     ######
###########################
Mymod <- "Mbin13A."
Saving(get(paste0(Mymod,num)))
num <- "26"
assign(paste0(Mymod,num),fitme(sp.mort.bin ~  sqrtBAIj.plot.1.mean + logdbh.plot.mean + 
                                      treeNbr + bio14_climate_mean.30 + 
                                      mean_spei12 + sqrtBA.ha.plot.1 +
                                      bio1_climate_mean.30:bio14_climate_mean.30 + 
                                      bio1_climate_mean.30:min_spei12 + bio1_climate_mean.30:mean_spei12 + 
                                      bio14_climate_mean.30:logBAj.plot.1 + 
                                      bio14_climate_mean.30:Plotcat +
                                      logBAj.plot.1:Plotcat + (1 | country), data=dfplot2,family=binomial,method='REML'))

MyMod <- get(paste0(Mymod,num))
MyMod <- as.data.frame(summary(MyMod)$beta_table)
MyMod2 <- MyMod[order(abs(MyMod$'t-value')),]
MyMod
MyMod2

###########################
#####   Bootstraps   ######
###########################
setwd(Dir) # Set dir in the new affined
### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[77]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[2], ignore.case = FALSE,fixed = T),get(load(file = list.files()[2]))) 
rm("x")
summary(get(sub(".rda","",list.files()[2])))

Diagnostic(Mbin13B.26,0.66,F)
Extraction(Mbin13B.26)
summary(Mbin13B.26)
ModelBoot(Mbin13B.26,4,10,LvL=20,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=10)
ModelBoot(Mbin13B.26,9,10,LvL=20,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=10)


#######################
##### AC Spatial ######  # Try to add spatial AC if possible
#######################
spaMM.options(separation_max=6363)
Mbin14A.19.AC <- fitme(sp.mort.bin ~ sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + logdbh.plot.mean + 
                         treeNbr + sqrtBA.ha.plot.1 + sqrtBA.O.plot.1 + Plotcat + 
                         I(bio14_climate_mean.30^2) + I(min_spei12^2) + tmean.djf_climate_mean.30:bio14_climate_mean.30 + 
                         tmean.djf_climate_mean.30:min_spei12 + tmean.djf_climate_mean.30:sqrtBA.ha.plot.1 + 
                         min_spei12:sqrtBA.O.plot.1 + mean_spei12:sqrtBA.ha.plot.1 + 
                         mean_spei12:sqrtBA.O.plot.1 + sqrtBA.ha.plot.1:sqrtBA.O.plot.1 + 
                         min_spei12:Plotcat + mean_spei12:Plotcat + (1 | country) + Matern(1|latitude + longitude),
                         data=dfplot2, family = binomial,method='REML')


###########################
##### ZT to affined  ###### # Same process with the ZT model 
###########################
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/Negbin"))
setwd(Dir)
### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[23]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[2], ignore.case = FALSE,fixed = T),get(load(file = list.files()[2]))) 
rm("x")
summary(get(sub(".rda","",list.files()[2])))
detach("package:mgcv", unload=TRUE)

Mymod <- "MnbZT13B."
Saving(get(paste0(Mymod,num)))
num <- "23"
assign(paste0(Mymod,num),fitme(sp.mortality.plot.rate.yr ~ sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + 
                                 min_spei12 + sqrtBA.ha.plot.1 + 
                                 logBAj.plot.1 + Plotcat +
                                 I(min_spei12^2) +
                                 bio1_climate_mean.30:sqrtBA.ha.plot.1 + bio1_climate_mean.30:logBAj.plot.1 + 
                                 bio14_climate_mean.30:min_spei12 + 
                                 mean_spei12:logBAj.plot.1 + 
                                 sqrtBA.ha.plot.1:logBAj.plot.1 +
                                 min_spei12:Plotcat + mean_spei12:Plotcat + 
                                 logBAj.plot.1:Plotcat + (1 | country),
                                 data=subset(dfplot2,sp.mortality.plot.rate.yr>0),family=negbin(),method='REML'))

MyMod <- get(paste0(Mymod,num))
MyMod <- as.data.frame(summary(MyMod)$beta_table)
MyMod2 <- MyMod[order(abs(MyMod$'t-value')),]
MyMod
#MyMod2


###########################
#####   Bootstraps   ######  ## ZT
###########################
setwd(Dir) # Set dir in the new affined
### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[75]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[2], ignore.case = FALSE,fixed = T),get(load(file = list.files()[2]))) 
rm("x")
summary(get(sub(".rda","",list.files()[2])))
Diagnostic(M2nbZT1A,0.66,F)
Extraction(MnbZT13B.23)
summary(MnbZT13B.23)
ModelBoot(MnbZT13B.23,5,6,LvL=20,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=10)
ModelBoot(MnbZT13B.23,3,6,LvL=20,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=10)
ModelBoot(MnbZT13B.23,9,6,LvL=20,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=10)

#######################
##### AC Spatial ######  # Try to add spatial AC if possible
#######################
MnbZT13B.23.AC <- fitme(sp.mortality.plot.rate.yr ~ sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + 
                          min_spei12 + sqrtBA.ha.plot.1 + 
                          logBAj.plot.1 + Plotcat +
                          I(min_spei12^2) +
                          bio1_climate_mean.30:sqrtBA.ha.plot.1 + bio1_climate_mean.30:logBAj.plot.1 + 
                          bio14_climate_mean.30:min_spei12 + 
                          mean_spei12:logBAj.plot.1 + 
                          sqrtBA.ha.plot.1:logBAj.plot.1 +
                          min_spei12:Plotcat + mean_spei12:Plotcat + 
                          logBAj.plot.1:Plotcat + (1 | country) + Matern(1|latitude + longitude),
                        data=subset(dfplot2,sp.mortality.plot.rate.yr>0),family=negbin(),method='REML')
Saving(MnbZT13B.23.AC)

