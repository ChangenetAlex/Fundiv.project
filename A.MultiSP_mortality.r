### Script sylvain 
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
library(dplyr)


Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
Allcode <- c("FAGSYL","PINHAL","PINSYL","QUEILE","ABIALB","QUEROB")
Allseuil <- c(0.7,0.7,0.8,0.8,0.7,0.8)
i = 1
# Here we kept only the variables we wanted 
for (i in 1:length(Allcode)){
    Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP")) # Directory 
    setwd(Dir)
    assign(paste0("dfplot",Allcode[i]),readRDS(paste0("dfplot",Allcode[i],Allseuil[i],".rds"))) #Base de données plot
    assign(paste0("dfplot2",Allcode[i]),get(paste0("dfplot",Allcode[i]))[,c("plotcode","latitude","longitude","speciesrichness","binomial","country","sp.mortality.plot.rate","mortality.plot.rate","sp.mortality.plot.rate.yr","sp.mortality.plot.count.yr","yearsbetweensurveys")]) #Base de données plot
    assign(paste0("dfplot3",Allcode[i]),get(paste0("dfplot2",Allcode[i]))[get(paste0("dfplot2",Allcode[i]))$sp.mortality.plot.rate.yr>0 & !is.na(get(paste0("dfplot2",Allcode[i]))$sp.mortality.plot.rate.yr),]) #Base de données plot
}

# For each species, transform continuous mortality and latitude as class variables made of three levels. 

dfplot3ABIALB$cut.latitude <- cut(dfplot3ABIALB$latitude, quantile(dfplot3ABIALB$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3ABIALB$cut.latitude) <- c(1,2,3)
dfplot3ABIALB$cut.sp.mortality.plot.rate.yr <- cut(dfplot3ABIALB$sp.mortality.plot.rate.yr, quantile(dfplot3ABIALB$sp.mortality.plot.rate.yr,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3ABIALB$cut.sp.mortality.plot.rate.yr) <- c(1,2,3)
summary(as.factor(dfplot3ABIALB$cut.sp.mortality.plot.rate.yr))
summary(as.factor(dfplot3ABIALB$cut.latitude))

dfplot3FAGSYL$cut.latitude <- cut(dfplot3FAGSYL$latitude, quantile(dfplot3FAGSYL$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3FAGSYL$cut.latitude) <- c(1,2,3)
dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr <- cut(dfplot3FAGSYL$sp.mortality.plot.rate.yr, quantile(dfplot3FAGSYL$sp.mortality.plot.rate.yr,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr) <- c(1,2,3)
summary(as.factor(dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr))
summary(as.factor(dfplot3FAGSYL$cut.latitude))

dfplot3PINHAL$cut.latitude <- cut(dfplot3PINHAL$latitude, quantile(dfplot3PINHAL$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3PINHAL$cut.latitude) <- c(1,2,3)
dfplot3PINHAL$cut.sp.mortality.plot.rate.yr <- cut(dfplot3PINHAL$sp.mortality.plot.rate.yr, quantile(dfplot3PINHAL$sp.mortality.plot.rate.yr,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3PINHAL$cut.sp.mortality.plot.rate.yr) <- c(1,2,3)
summary(as.factor(dfplot3PINHAL$cut.sp.mortality.plot.rate.yr))
summary(as.factor(dfplot3PINHAL$cut.latitude))

dfplot3PINSYL$cut.latitude <- cut(dfplot3PINSYL$latitude, quantile(dfplot3PINSYL$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3PINSYL$cut.latitude) <- c(1,2,3)
dfplot3PINSYL$cut.sp.mortality.plot.rate.yr <- cut(dfplot3PINSYL$sp.mortality.plot.rate.yr, quantile(dfplot3PINSYL$sp.mortality.plot.rate.yr,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3PINSYL$cut.sp.mortality.plot.rate.yr) <- c(1,2,3)
summary(as.factor(dfplot3PINSYL$cut.sp.mortality.plot.rate.yr))
summary(as.factor(dfplot3PINSYL$cut.latitude))

dfplot3QUEROB$cut.latitude <- cut(dfplot3QUEROB$latitude, quantile(dfplot3QUEROB$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3QUEROB$cut.latitude) <- c(1,2,3)
dfplot3QUEROB$cut.sp.mortality.plot.rate.yr <- cut(dfplot3QUEROB$sp.mortality.plot.rate.yr, quantile(dfplot3QUEROB$sp.mortality.plot.rate.yr,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3QUEROB$cut.sp.mortality.plot.rate.yr) <- c(1,2,3)
summary(as.factor(dfplot3QUEROB$cut.sp.mortality.plot.rate.yr))
summary(as.factor(dfplot3QUEROB$cut.latitude))

dfplot3QUEILE$cut.latitude <- cut(dfplot3QUEILE$latitude, quantile(dfplot3QUEILE$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3QUEILE$cut.latitude) <- c(1,2,3)
dfplot3QUEILE$cut.sp.mortality.plot.rate.yr <- cut(dfplot3QUEILE$sp.mortality.plot.rate.yr, quantile(dfplot3QUEILE$sp.mortality.plot.rate.yr,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3QUEILE$cut.sp.mortality.plot.rate.yr) <- c(1,2,3)
summary(as.factor(dfplot3QUEILE$cut.sp.mortality.plot.rate.yr))
summary(as.factor(dfplot3QUEILE$cut.latitude))


## For each species, sample randomly three plots for each combination of modality. 
# (3 levels for latitude, 3 levels for mortality => 9 modality * 3 plots = 27 rows.)

table(dfplot3ABIALB$cut.latitude,dfplot3ABIALB$cut.sp.mortality.plot.rate.yr)
k=1
for (i in 1:3){
  for (j in 1:3){
    assign(paste0("test",k),sample_n(dfplot3ABIALB[dfplot3ABIALB$cut.latitude==i & dfplot3ABIALB$cut.sp.mortality.plot.rate.yr==j,],3))
    k = k+1
    j=j+1
    }
  i=i+1
}
dfplot4ABIALB <- rbind(test1,test2,test3,test4,test5,test6,test7,test8,test9)

table(dfplot3FAGSYL$cut.latitude,dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr)
k=1
for (i in 1:3){
  for (j in 1:3){
    assign(paste0("test",k),sample_n(dfplot3FAGSYL[dfplot3FAGSYL$cut.latitude==i & dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr==j,],3))
    k = k+1
    j=j+1
  }
  i=i+1
}
dfplot4FAGSYL <- rbind(test1,test2,test3,test4,test5,test6,test7,test8,test9)


table(dfplot3PINHAL$cut.latitude,dfplot3PINHAL$cut.sp.mortality.plot.rate.yr)
k=1
for (i in 1:3){
  for (j in 1:3){
    assign(paste0("test",k),sample_n(dfplot3PINHAL[dfplot3PINHAL$cut.latitude==i & dfplot3PINHAL$cut.sp.mortality.plot.rate.yr==j,],3))
    k = k+1
    j=j+1
  }
  i=i+1
}
dfplot4PINHAL <- rbind(test1,test2,test3,test4,test5,test6,test7,test8,test9)


table(dfplot3PINSYL$cut.latitude,dfplot3PINSYL$cut.sp.mortality.plot.rate.yr)
k=1
for (i in 1:3){
  for (j in 1:3){
    assign(paste0("test",k),sample_n(dfplot3PINSYL[dfplot3PINSYL$cut.latitude==i & dfplot3PINSYL$cut.sp.mortality.plot.rate.yr==j,],3))
    k = k+1
    j=j+1
  }
  i=i+1
}
dfplot4PINSYL <- rbind(test1,test2,test3,test4,test5,test6,test7,test8,test9)


table(dfplot3QUEROB$cut.latitude,dfplot3QUEROB$cut.sp.mortality.plot.rate.yr)
k=1
for (i in 1:3){
  for (j in 1:3){
    assign(paste0("test",k),sample_n(dfplot3QUEROB[dfplot3QUEROB$cut.latitude==i & dfplot3QUEROB$cut.sp.mortality.plot.rate.yr==j,],3))
    k = k+1
    j=j+1
  }
  i=i+1
}
dfplot4QUEROB <- rbind(test1,test2,test3,test4,test5,test6,test7,test8,test9)


table(dfplot3QUEILE$cut.latitude,dfplot3QUEILE$cut.sp.mortality.plot.rate.yr)
k=1
for (i in 1:3){
  for (j in 1:3){
    assign(paste0("test",k),sample_n(dfplot3QUEILE[dfplot3QUEILE$cut.latitude==i & dfplot3QUEILE$cut.sp.mortality.plot.rate.yr==j,],3))
    k = k+1
    j=j+1
  }
  i=i+1
}
dfplot4QUEILE <- rbind(test1,test2,test3,test4,test5,test6,test7,test8,test9)



saveRDS(dfplot4ABIALB, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4ABIALB.rds"))
saveRDS(dfplot4FAGSYL, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4FAGSYL.rds"))
saveRDS(dfplot4PINHAL, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4PINHAL.rds"))
saveRDS(dfplot4PINSYL, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4PINSYL.rds"))
saveRDS(dfplot4QUEILE, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4QUEILE.rds"))
saveRDS(dfplot4QUEROB, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4QUEROB.rds"))

####################################
### Partie 2 multispecies sites ! ##
####################################

## 3 species 
library(plyr)
elem <- c("dfplot3ABIALB","dfplot3FAGSYL","dfplot3PINHAL","dfplot3PINSYL","dfplot3QUEILE","dfplot3QUEROB")
t(combn(elem, 3))
# These lines are all combinations for which plots containing all three species were found. 
df3sp1 <- join_all(list(dfplot3ABIALB,dfplot3FAGSYL,dfplot3PINSYL), by='plotcode', type='inner') #2
df3sp1[,c(9,21,33)]
df3sp2 <- join_all(list(dfplot3ABIALB,dfplot3PINSYL,dfplot3QUEROB), by='plotcode', type='inner') #1
df3sp2[,c(9,21,33)]
df3sp3 <- join_all(list(dfplot3FAGSYL,dfplot3PINSYL,dfplot3QUEROB), by='plotcode', type='inner') #2
df3sp3[,c(9,21,33)]
df3sp4 <- join_all(list(dfplot3PINHAL,dfplot3PINSYL,dfplot3QUEILE), by='plotcode', type='inner') #3
df3sp4[,c(9,21,33)]
df3sp5 <- join_all(list(dfplot3PINHAL,dfplot3QUEILE,dfplot3QUEROB), by='plotcode', type='inner') #1
df3sp5[,c(9,21,33)]

df3sp <- rbind(df3sp1,df3sp2,df3sp3,df3sp4,df3sp5)
saveRDS(df3sp, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/df3sp.rds"))


## 2 species
## Check all combinations
## Create a new variable dif.mort that is the absolute difference in mortality between the two species for each plot. 
## Among all possible plots : select the 10 for which this values is the highest (when thre are more than 10 plots)

df2sp1 <- join_all(list(dfplot3ABIALB,dfplot3FAGSYL), by='plotcode', type='inner') #68 plots 
df2sp1$dif.mort <- abs(df2sp1[,c(9)]-df2sp1[,c(21)])
df2sp1 <- df2sp1[order(-df2sp1$dif.mort)[1:10],]

df2sp2 <- join_all(list(dfplot3ABIALB,dfplot3PINSYL), by='plotcode', type='inner') #58
df2sp2$dif.mort <- abs(df2sp2[,c(9)]-df2sp2[,c(21)])
df2sp2 <- df2sp2[order(-df2sp2$dif.mort)[1:10],]

df2sp3 <- join_all(list(dfplot3ABIALB,dfplot3QUEROB), by='plotcode', type='inner') #8
df2sp3$dif.mort <- abs(df2sp3[,c(9)]-df2sp3[,c(21)])

df2sp4 <- join_all(list(dfplot3FAGSYL,dfplot3PINSYL), by='plotcode', type='inner') #66
df2sp4$dif.mort <- abs(df2sp4[,c(9)]-df2sp4[,c(21)])
df2sp4 <- df2sp4[order(-df2sp4$dif.mort)[1:10],]

df2sp5 <- join_all(list(dfplot3FAGSYL,dfplot3QUEROB), by='plotcode', type='inner') #31
df2sp5$dif.mort <- abs(df2sp5[,c(9)]-df2sp5[,c(21)])
df2sp5 <- df2sp5[order(-df2sp5$dif.mort)[1:10],]

df2sp6 <- join_all(list(dfplot3PINHAL,dfplot3PINSYL), by='plotcode', type='inner') #21
df2sp6$dif.mort <- abs(df2sp6[,c(9)]-df2sp6[,c(21)])
df2sp6 <- df2sp6[order(-df2sp6$dif.mort)[1:10],]

df2sp7 <- join_all(list(dfplot3PINHAL,dfplot3QUEILE), by='plotcode', type='inner') #47
df2sp7$dif.mort <- abs(df2sp7[,c(9)]-df2sp7[,c(21)])
df2sp7 <- df2sp7[order(-df2sp7$dif.mort)[1:10],]

df2sp8 <- join_all(list(dfplot3PINHAL,dfplot3QUEROB), by='plotcode', type='inner') #1
df2sp8$dif.mort <- abs(df2sp8[,c(9)]-df2sp8[,c(21)])

df2sp9 <- join_all(list(dfplot3PINSYL,dfplot3QUEROB), by='plotcode', type='inner') #62
df2sp9$dif.mort <- abs(df2sp9[,c(9)]-df2sp9[,c(21)])
df2sp9 <- df2sp9[order(-df2sp9$dif.mort)[1:10],]

df2sp10 <- join_all(list(dfplot3PINSYL,dfplot3QUEILE), by='plotcode', type='inner') #27
df2sp10$dif.mort <- abs(df2sp10[,c(9)]-df2sp10[,c(21)])
df2sp10 <- df2sp10[order(-df2sp10$dif.mort)[1:10],]

df2sp11<- join_all(list(dfplot3QUEILE,dfplot3QUEROB), by='plotcode', type='inner') #4
df2sp11$dif.mort <- abs(df2sp11[,c(9)]-df2sp11[,c(21)])


df2sp <- rbind(df2sp1,df2sp2,df2sp3,df2sp4,df2sp5,df2sp6,df2sp7,df2sp8,df2sp9,df2sp10,df2sp11)
saveRDS(df2sp, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/df2sp.rds"))


    