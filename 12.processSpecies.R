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

###################################################################################
###                                                                            ####
###                                                                            #### # Here, 7 scripts are run to obtain the finals global databases
###     PART one : how to extract the information we want for all species      ####
###                                                                            ####
###                                                                            ####
###################################################################################




############################################################
####    Extraction de la base de donnée climatique     #####
############################################################

Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/6.climateSpecies.R"))
Allcode <- c("PINSYL")
#Allcode <- c("PICABI","PINPINA","FAGSYL","PINHAL","QUEROB","PINNIG","QUEPET","CASSAT","ABIALB","QUEPUB","QUEPYR","FRAEXC","PINPIN",
             #"QUESUB","PINRAD","BETPEN","PINCAN","PINUNC","EUCGLO","QUEFAG","ALNGLU","POPTRE","ACEPSE","LARDEC","ROBPSE","POPNIG","PICSIT","QUERUB")
for (code in Allcode){
  try(ExtractClimate(CODE=code,yearINI=yearI,yearFIN=yearF,Intervalle=Inter,files.wanted=files.w,save=T),silent=T) # Maybe not all the intervalle
}

############################################################
####       Spatialisation of the wanted species        #####
############################################################
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/8.spatialisation.R"))
Allcode <- c("PINSYL")
#Allcode <- c("PINHAL","QUEROB","PINNIG","QUEPET","CASSAT","ABIALB","QUEPYR","FRAEXC","PINPIN",
#             "QUESUB","BETPEN")
#Allcode <- c("FAGSYL","PINHAL","QUEROB","PINNIG","QUEPET","CASSAT","ABIALB","QUEPUB","QUEPYR","FRAEXC","PINPIN",
     #        "QUESUB","BETPEN","ALNGLU","POPTRE","ACEPSE","LARDEC","POPNIG")
Errors.files <- file(paste0(Dir,"species/Errors.spatialisation.Rout"), open="wt")
sink(file=Errors.files,append = T, type="message")
for (code in Allcode){
  try(Spatialisation(CODE = code,SaveAll = T,myPointsAll = T),silent=F) #re run the code with this arguments 
}
sink(type="message")
close(Errors.files)



############################################################
####    Cut marginality for all species and do it again  ###
####        for species whose plot number is too low     ###
############################################################

Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/9.Marginality.R"))
Allcode <- c("PINSYL")
#Allcode <- list("PICABI","PINPINA","QUEILE","POPTRE","LARDEC","POPNIG")
Errors.files2 <- file(paste0(Dir,"species/Errors.Marginality.Rout"), open="wt")
sink(Errors.files2, type="message")
for (code in Allcode){
  try(Acp1000(CODE = code,nsample=10000,NF=2,seuil=0.8,seuilC=0.6,save=T),silent=T) #re run the code with this arguments # First function for PCA
}
#Allcode <- c("PICABI","PINPINA","FAGSYL","PINHAL","QUEROB","QUEILE","PINNIG","QUEPET","CASSAT","ABIALB","QUEPUB","QUEPYR","FRAEXC","PINPIN",
#             "QUESUB","BETPEN","ALNGLU","POPTRE","ACEPSE","LARDEC","POPNIG")
for (code in Allcode){
  try(Marginality_Levels(CODE=code,seuil=0.8,seuilC=0.6,save=T),silent=T)
}
sink(type="message")
close(Errors.files2)




############################################################
####    SPEI indexes calculation for all Intervalles   #####
############################################################

Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/10bis.SPEI.R"))
Years <- c("spei01","spei03","spei06","spei18","spei24","spei36","spei48")
for (years in Years){
  try(SPEI(nStart = nStart,nEnd = nEnd, ncname = years),silent=T)
}


############################################################
####    Mean DBH and mean BAIj.plot & cumulated DBH    #####
############################################################

Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/10bis.meanDBH.R"))
Allcode <- c("PINSYL","FAGSYL","PICABI","PINPINA","PINHAL","QUEROB","QUEILE","PINNIG","QUEPET","CASSAT","ABIALB","QUEPUB","QUEPYR","FRAEXC","PINPIN",
             "QUESUB","BETPEN","ALNGLU","POPTRE","ACEPSE","LARDEC","POPNIG")
Errors.files <- file(paste0(Dir,"our-data/species/Errors.dfplots.meanDBH.Rout"), open="wt")
sink(Errors.files, type="message")
for (code in Allcode){
  for (s in c(0.8,0.7)){
  try(MeanDBH(dir="bureau",CODE = code,seuil=s,seuilC=0.6),silent=T) #re run the code with this arguments # First function for PCA
}}
sink(type="message")
close(Errors.files)



############################################################
####        From individuals trees to plot data        #####
############################################################

Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/10.Dfplots.R"))
Allcode <- list("FAGSYL")
Allcode <- c("PINSYL","FAGSYL","PICABI","PINPINA","PINHAL","QUEROB","QUEILE","PINNIG","QUEPET","CASSAT","ABIALB","QUEPUB","QUEPYR","FRAEXC","PINPIN",
             "QUESUB","BETPEN","ALNGLU","POPTRE","ACEPSE","LARDEC","POPNIG")
#Allcode <- "PINSYL"
Errors.files <- file(paste0(Dir,"our-data/species/Errors.dfplots.Rout"), open="wt")
sink(Errors.files,type="message")
for (code in Allcode){
  for (s in c(0.8,0.7)){
  try(MyDFs(dir="bureau",CODE = code,seuil=s,seuilC=0.6),silent=T) #re run the code with this arguments # First function for PCA
}}
sink(type="message")
close(Errors.files)

CODE = "FAGSYL"
seuil = 0.7

Saving(Mbin1)
