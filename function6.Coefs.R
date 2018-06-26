Mbin1FAGSYLspaMM <- fitme(sp.mort.bin ~ 0 + treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1 + bio14_climate_mean.30 + logbio14 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1 + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:logbio14 + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:logbio14 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1:logbio14 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:logbio1 + Plotcat:logbio14 + (1|country), data=dfplot2,family = binomial(),method='REML')
Mbin1FAGSYLspaMM <- fitme(sp.mort.bin ~ treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1_climate_mean.30 + bio14_climate_mean.30 + logbio14_climate_mean.30 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1_climate_mean.30 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1_climate_mean.30 + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:logbio14_climate_mean.30 + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:logbio14_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1_climate_mean.30:logbio14_climate_mean.30 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:logbio1_climate_mean.30 + Plotcat:logbio14_climate_mean.30 + (1|country), data=dfplot2,family = binomial(),method='REML')
Saving(Mbin1FAGSYLspaMM)


# Dans mon modèle, les paramètres d'interactions doivent être après les para log, seuls et cubique
# Les noms des variables scaled doivent être les même que les non scales et pas une abbréviation

# get the parameter values from the glm model Unique parameters of the model
library(reshape2)


ExtracTest <- function(x){
        A <- paste0(getCall(x)[2])
        A <- unlist(strsplit(A, "~",fixed = T))[2]
        A <- unlist(strsplit(A, "+",fixed = T))
        A <- grep(A,pattern = "|",fixed=T,value=T,invert=T)
        A <- grep(A,pattern = ":",fixed=T,value=T,invert=T)
        A <- grep(A,pattern = "e)",fixed=T,value=T,invert=T) #Remove the AC term which remains
        A <- sub("I(","", A, ignore.case = FALSE,fixed = T)
        A <- sub("^2)","", A, ignore.case = FALSE,fixed = T)
        A <- sub("offset(log(","", A, ignore.case = FALSE,fixed = T)
        A <- sub("))","", A, ignore.case = FALSE,fixed = T)
        A <- sub("\n ","", A, ignore.case = FALSE,fixed = T) # New line
        A <- sub("  ","", A, ignore.case = FALSE,fixed = T)
        A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
        A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
        assign("A",unique(A),envir = .GlobalEnv) 
        A <- unique(A)
        A
}
Effect_coef <- function(x,y){
  if (grepl(deparse(substitute(x)),pattern="bin",fixed=T)==T){
    Dir <- paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/FAGSYL/CLIMAP/Models/binomial/",deparse(substitute(x)),"/")
  }else Dir <- paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/FAGSYL/CLIMAP/Models/Negbin/",deparse(substitute(x)),"/")
  setwd(Dir)      
  s <- summary(x)
        if(length(list.files(path=Dir,pattern=paste0("clim_effect_",deparse(substitute(x)),".RData")))==0){
                clim.effect <- dfplot2[,c(A,"latitude")]
        }else {clim.effect <- get(load(paste0("clim_effect_",deparse(substitute(x)),".RData")))}
        # For any parameters 
        B <- rownames(s$beta_table)
        B <- grep(B,pattern = y,fixed=T,value=T,invert=F)
        B <- grep(B,pattern = "Plotcat",fixed=T,value=T,invert=T) # To check remove the plotcat
        # Here B is the list of the coefs names
        C1 <- grep(B,pattern = ":",fixed=T,value=T,invert=T) # No interactions effects
        C2 <- grep(B,pattern = ":",fixed=T,value=T,invert=F) # Interactions effects
        C2 <- unlist(strsplit(C2, ":",fixed = T))
        C2 <- grep(C2,pattern = y,fixed=T,value=T,invert=T) # interactions effects
        
        clim.effect[,paste0("effect_",y)] <- sum(s$beta_table[B[1:length(C1)],1]) # Ici on fait l'hypothèse que les coefs simple sont dans l'ordre avant les interactions. Je dois me débrouiller pour que ce soit vrai tout le temps  
        if (length(C2)>0){ # If just single effect or interactions 
        for (i in 1:length(C2)){
                clim.effect[,paste0("effect_",y)] <- clim.effect[,paste0("effect_",y)] + s$beta_table[B[i+length(C1)],1]*clim.effect[,C2[i]]
        }}
        clim.effect <- as.data.frame(clim.effect)
        save(clim.effect, file=paste0("clim_effect_",deparse(substitute(x)),".RData"))
        print(clim.effect[1:10,]) # Check if I have all the columns I want 
        }
     
ExtracTest(Mbin1FAGSYLspaMM) # Me donne la liste des paramètres de mon modèles. 
for (i in c("BAIj.plot.bis","BA.ha.plot.1","BAj.plot.1","treeNbr","bio1_climate_mean.30","bio14_climate_mean.30")){
  Effect_coef(Mbin1FAGSYLspaMM,i)}     # Para que l'on veut regarder   



### 3eme function to write 

### Load my Df ###
clim.effect <- get(load(paste0("clim_effect_",deparse(substitute(Mbin1FAGSYLspaMM)),".RData")))
EffectCol <- grep(colnames(clim.effect),pattern = "effect",fixed=T,value=T,invert=F) # Extract the column for which the coef were extracted
EffectCum <- c("climate","biotic","competition") # Add the three categories 
EffectAll <- c(EffectCol,EffectCum)
# Newdataframe made of these extracted effect and the latitude
ind.clim.bio.abs.imp <- as.data.frame(clim.effect[,c("latitude")])
names(ind.clim.bio.abs.imp) <- c("latitude")

ind.clim.bio.rel.imp <- as.data.frame(clim.effect[,c("latitude")])
names(ind.clim.bio.rel.imp) <- c("latitude")

# As many varaibles as the number of runs of the previous function and fill it with NA
for (i in 1:length(EffectCol)){
ind.clim.bio.rel.imp[,EffectCol[i]] <- NA
ind.clim.bio.abs.imp[,EffectCol[i]] <- NA
}
# Extract the max OR the sum for each latitude
Var <- grep(colnames(clim.effect),pattern = "effect",fixed=T,value=T,invert=F) # Plus de log et termes au carré
Var.comp <- grep(Var,pattern = paste(c("BA","tree"),collapse="|"),value=T,invert=F)
Var.comp <- grep(Var.comp,pattern = "BAI",fixed=T,value=T,invert=T) # Variable de compétition tree numbers as well
Var.climate <- grep(Var,pattern = paste(c("spei","climate"),collapse="|"),value=T,invert=F) # Variable de climat (spei + worlclim)
Var.biotic <- grep(Var,pattern = paste(c("BA","tree"),collapse="|"),value=T,invert=F) # Variable compétition + BAI

if (length(Var.comp)>1){clim.effect$competition <- apply(abs(clim.effect[,Var.comp]), 1, function(x) sum(x)/length(Var.comp)) #need length superior to 1
}else clim.effect$competition <- abs(clim.effect[,Var.comp])

if (length(Var.climate)>1){clim.effect$climate <- apply(abs(clim.effect[,Var.climate]), 1, function(x) sum(x)/length(Var.climate)) #need length superior to 1
}else clim.effect$climate <- abs(clim.effect[,Var.climate])

if (length(Var.biotic)>1){clim.effect$biotic <- apply(abs(clim.effect[,Var.biotic]), 1, function(x) sum(x)/length(Var.biotic)) #need length superior to 1
}else clim.effect$biotic <- abs(clim.effect[,Var.biotic])

clim.effect$max <- apply(cbind(abs(clim.effect[,EffectCol])), 1, max) # Max 
clim.effect$sum <- apply(cbind(abs(clim.effect[,EffectCol])), 1, sum) # Sum 
clim.effect$maxcum <- apply(cbind(abs(clim.effect[,EffectCum])), 1, max) # Max 
clim.effect$sumcum <- apply(cbind(abs(clim.effect[,EffectCum])), 1, sum) # Sum 

for (i in 1:length(EffectCol)){
  ind.clim.bio.abs.imp[,EffectCol[i]] <- abs(clim.effect[,EffectCol[i]])/clim.effect$max  # Importance compare to the max one
  ind.clim.bio.rel.imp[,EffectCol[i]] <- abs(clim.effect[,EffectCol[i]])/clim.effect$sum  # Relative importance compare to the others (sum is qual to one)
}
for (i in 1:length(EffectCum)){
  ind.clim.bio.abs.imp[,EffectCum[i]] <- abs(clim.effect[,EffectCum[i]])/clim.effect$maxcum  # Importance compare to the max one
  ind.clim.bio.rel.imp[,EffectCum[i]] <- abs(clim.effect[,EffectCum[i]])/clim.effect$sumcum  # Relative importance compare to the others (sum is qual to one)
}


save(ind.clim.bio.abs.imp, file="./ind_abs_imp_all_variables.RData") # to modify 
save(ind.clim.bio.rel.imp, file="./ind_rel_imp_all_variables.RData") # to modify
#load("./ind_rel_imp_all_variables.RData")


ind.clim.bio.rel.imp$latitude_c <- cut(ind.clim.bio.rel.imp$latitude, seq(floor(min(ind.clim.bio.rel.imp$latitude)),ceiling(max(ind.clim.bio.rel.imp$latitude)),0.5))
ind.clim.bio.abs.imp$latitude_c <- cut(ind.clim.bio.abs.imp$latitude, seq(floor(min(ind.clim.bio.abs.imp$latitude)),ceiling(max(ind.clim.bio.abs.imp$latitude)),0.5))

t.mean.rel <- ddply(ind.clim.bio.rel.imp, .(latitude_c), colwise(mean,c(EffectAll)), .drop=F,na.rm=T) #ok 
t.max.rel <- ddply(ind.clim.bio.rel.imp, .(latitude_c), colwise(max,na.rm=T,c(EffectAll)),.drop=F,na.rm=T,.inform = T) # Pbm with NA
t.sd.rel <- ddply(ind.clim.bio.rel.imp, .(latitude_c), colwise(sd,c(EffectAll)), .drop=F,na.rm=T) #ok
t.n.rel <- ddply(ind.clim.bio.rel.imp, .(latitude_c), colwise(length,c(EffectAll)), .drop=F) #ok

t.mean.abs <- ddply(ind.clim.bio.abs.imp, .(latitude_c), colwise(mean,c(EffectAll)), .drop=F,na.rm=T) #ok 
t.max.abs <- ddply(ind.clim.bio.abs.imp, .(latitude_c), colwise(max,na.rm=T,c(EffectAll)),.drop=F,na.rm=T,.inform = T) # Pbm with NA
t.sd.abs <- ddply(ind.clim.bio.abs.imp, .(latitude_c), colwise(sd,c(EffectAll)), .drop=F,na.rm=T) #ok
t.n.abs <- ddply(ind.clim.bio.abs.imp, .(latitude_c), colwise(length,c(EffectAll)), .drop=F) #ok

for (i in 1:length(EffectAll)){
  t.sd.rel[,paste0(EffectAll[i],"_se")] <- t.sd.rel[,EffectAll[i]]/sqrt(t.n.rel[,EffectAll[i]]) # Calcul of standard errors 
  t.sd.abs[,paste0(EffectAll[i],"_se")] <- t.sd.abs[,EffectAll[i]]/sqrt(t.n.abs[,EffectAll[i]]) # Calcul of standard errors 
}

lats <- seq(floor(min(ind.clim.bio.rel.imp$latitude)),ceiling(max(ind.clim.bio.rel.imp$latitude)),0.5)[-1] # Remove the first 
t.mean.rel$latitude <- lats
t.sd.rel$latitude <- lats
t.mean.abs$latitude <- lats
t.sd.abs$latitude <- lats

ind.rel.imp.clim.bio.long <- melt(t.mean.rel, id = "latitude", measure = c(EffectAll))
names(ind.rel.imp.clim.bio.long) <- c("latitude", "variable", "mean")
ind.abs.imp.clim.bio.long <- melt(t.mean.abs, id = "latitude", measure = c(EffectAll))
names(ind.abs.imp.clim.bio.long) <- c("latitude", "variable", "mean")
## error bands
t.se.rel <- t.sd.rel[,(ncol(t.sd.rel)-length(EffectAll)):ncol(t.sd.rel)]
colnames(t.se.rel) <- sub("_se","", colnames(t.se.rel), ignore.case = FALSE,fixed = T)
ind.rel.imp.clim.bio.long.se <- melt(t.se.rel, id = "latitude", measure = c(EffectAll))
names(ind.rel.imp.clim.bio.long.se) <- c("latitude", "variable", "se")

t.se.abs <- t.sd.abs[,(ncol(t.sd.abs)-length(EffectAll)):ncol(t.sd.abs)]
colnames(t.se.abs) <- sub("_se","", colnames(t.se.abs), ignore.case = FALSE,fixed = T)
ind.abs.imp.clim.bio.long.se <- melt(t.se.abs, id = "latitude", measure = c(EffectAll))
names(ind.abs.imp.clim.bio.long.se) <- c("latitude", "variable", "se")


# merge the mean and the se dataframes
clim.bio.rel <- merge(ind.rel.imp.clim.bio.long, ind.rel.imp.clim.bio.long.se, by=c("latitude", "variable"))
clim.bio.rel$lwr <- clim.bio.rel$mean-(1.96*clim.bio.rel$se)
clim.bio.rel$upr <- clim.bio.rel$mean+(1.96*clim.bio.rel$se) 
clim.bio.rel$lwr <- ifelse(clim.bio.rel$lwr<0,0,clim.bio.rel$lwr)
clim.bio.rel$upr <- ifelse(clim.bio.rel$upr>1,1,clim.bio.rel$upr)

clim.bio.abs <- merge(ind.abs.imp.clim.bio.long, ind.abs.imp.clim.bio.long.se, by=c("latitude", "variable"))
clim.bio.abs$lwr <- clim.bio.abs$mean-(1.96*clim.bio.abs$se)
clim.bio.abs$upr <- clim.bio.abs$mean+(1.96*clim.bio.abs$se) 
clim.bio.abs$lwr <- ifelse(clim.bio.abs$lwr<0,0,clim.bio.abs$lwr)
clim.bio.abs$upr <- ifelse(clim.bio.abs$upr>1,1,clim.bio.abs$upr)

save(clim.bio.rel, file=paste0("clim_bio_rel_",deparse(substitute(x)),".RData"))
save(clim.bio.abs, file=paste0("clim_bio_abs_",deparse(substitute(x)),".RData"))




### ggplot part = 4eme function  