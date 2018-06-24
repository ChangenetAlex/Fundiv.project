Mbin1FAGSYLspaMM <- fitme(sp.mort.bin ~ 0 + treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1 + bio14_climate_mean.30 + logbio14 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1 + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:logbio14 + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:logbio14 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1:logbio14 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:logbio1 + Plotcat:logbio14 + (1|country), data=dfplot2,family = binomial(),method='REML')
Mbin1FAGSYLspaMM <- fitme(sp.mort.bin ~ treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1 + bio14_climate_mean.30 + logbio14 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1 + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:logbio14 + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:logbio14 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1:logbio14 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:logbio1 + Plotcat:logbio14 + (1|country), data=dfplot2,family = binomial(),method='REML')


# Dans mon modèle, les paramètres d'interactions doivent être après les para log, seuls et cubique
# Les noms des variables scaled doivent être les même que les non scales et pas une abbréviation

# get the parameter values from the glm model

##### CALCULATE CLIMATIC AND COMPETITION EFFECT - From individual data ##############

# Store the individual-level data used in the model (standardised values are used in the model, 
# so standardised values are used when calculating the effect sizes)
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
ExtracTest(Mbin1FAGSYLspaMM) # Me donne la liste des paramètres de mon modèles. 

Effect_coef <- function(x,y){
        s <- summary(x)
        if(length(list.files(path=Dir,pattern="clim_effect.RData"))==0){
                clim.effect <- dfplot2[,c(A,"latitude")]
        }else {load("clim_effect.RData")}
        # For any parameters 
        B <- rownames(s$beta_table)
        B <- grep(B,pattern = y,fixed=T,value=T,invert=F)
        B <- grep(B,pattern = "Plotcat",fixed=T,value=T,invert=T) # To check remove the plotcat
        # Here B is the list of the coefs names
        C1 <- grep(B,pattern = ":",fixed=T,value=T,invert=T) # No interactions effects
        C2 <- grep(B,pattern = ":",fixed=T,value=T,invert=F) # Interactions effects
        C2 <- unlist(strsplit(C2, ":",fixed = T))
        C2 <- grep(C2,pattern = y,fixed=T,value=T,invert=T) # interactions effects
        
        clim.effect[,paste0("effect",y)] <- sum(s$beta_table[B[1:length(C1)],1]) # Ici on fait l'hypothèse que les coefs simple sont dans l'ordre avant les interactions. Je dois me débrouiller pour que ce soit vrai tout le temps  
        for (i in 1:length(C2)){
                clim.effect[,paste0("effect",y)] <- clim.effect[,paste0("effect",y)] + s$beta_table[B[i+length(C1)],1]*clim.effect[,C2[i]]
        }
        clim.effect <- as.data.frame(clim.effect)
        save(clim.effect, file="./clim_effect.RData")
        print(clim.effect[1:10,]) # Check if I have all the columns I want 
        }
     
ExtracTest(Mbin1FAGSYLspaMM) #Liste des para 
Effect_coef(Mbin1FAGSYLspaMM,"bio14")     # Para que l'on veut regarder   
load("clim_effect.RData")



library(reshape2)

setwd("/home/juliette/Documents/STAGES/Stage2018/a/mydata/pinus/data/")

##### 1) Relative importance of minSPEI, FROSTS, PET and competition - from individual data #################################

load("./clim_effect.RData")

ind.clim.bio.rel.imp <- as.data.frame(clim.effect[,c("latitude")])
names(ind.clim.bio.rel.imp) <- c("latitude")
ind.clim.bio.rel.imp$BAj.plot.1 <- NA
ind.clim.bio.rel.imp$bio14 <- NA

clim.effect$max <- apply(cbind(abs(clim.effect$effectBAj.plot.1), abs(clim.effect$effectbio14)), 1, max)

ind.clim.bio.rel.imp$BAj.plot.1 <- abs(clim.effect$effectBAj.plot.1)/clim.effect$max 
ind.clim.bio.rel.imp$bio14 <- abs(clim.effect$effectbio14)/clim.effect$max 

names(ind.clim.bio.rel.imp) <- c("Latitude", "minSPEI", "PET", "FROSTS","Competition")

save(ind.clim.bio.rel.imp, file="./ind_rel_imp_all_variables.RData")
load("./ind_rel_imp_all_variables.RData")

ind.clim.bio.rel.imp$Latitude_c <- cut(ind.clim.bio.rel.imp$Latitude, seq(36.5,70,0.5))

t.mean <- aggregate(ind.clim.bio.rel.imp[,c("minSPEI", "PET", "FROSTS","Competition")], by=list(ind.clim.bio.rel.imp$Latitude_c), mean)
t.max <- aggregate(ind.clim.bio.rel.imp[,c("minSPEI", "PET", "FROSTS","Competition")], by=list(ind.clim.bio.rel.imp$Latitude_c), max)
t.sd <- aggregate(ind.clim.bio.rel.imp[,c("minSPEI", "PET", "FROSTS","Competition")], by=list(ind.clim.bio.rel.imp$Latitude_c), sd)
t.n <- aggregate(ind.clim.bio.rel.imp[,c("minSPEI", "PET", "FROSTS","Competition")], by=list(ind.clim.bio.rel.imp$Latitude_c), length)

t.sd$minSPEI.se <- t.sd$minSPEI/sqrt(t.n$minSPEI)
t.sd$PET.se <- t.sd$PET/sqrt(t.n$PET)
t.sd$FROSTS.se <- t.sd$FROSTS/sqrt(t.n$FROSTS)
t.sd$Competition.se <- t.sd$Competition/sqrt(t.n$Competition)

lats <- seq(36.5,69.5,0.5)
lats <- lats[!lats %in% c(37.5,38,38.5,39)]
t.mean$Latitude <- lats
t.sd$Latitude <- lats

ind.rel.imp.clim.bio.long <- melt(t.mean, id = "Latitude", measure = c("minSPEI", "PET", "FROSTS","Competition"))
names(ind.rel.imp.clim.bio.long) <- c("Latitude", "variable", "mean")

## error bands
t.se <- t.sd[,6:10]
names(t.se) <- c("minSPEI", "PET", "FROSTS","Competition", "Latitude")
ind.rel.imp.clim.bio.long.se <- melt(t.se, id = "Latitude", measure = c("minSPEI", "PET", "FROSTS","Competition"))
names(ind.rel.imp.clim.bio.long.se) <- c("Latitude", "variable", "se")

# merge the mean and the se dataframes
clim.bio <- merge(ind.rel.imp.clim.bio.long, ind.rel.imp.clim.bio.long.se, by=c("Latitude", "variable"))
clim.bio$lwr <- clim.bio$mean-(1.96*clim.bio$se)
clim.bio$upr <- clim.bio$mean+(1.96*clim.bio$se) 

clim.bio$lwr <- ifelse(clim.bio$lwr<0,0,clim.bio$lwr)
clim.bio$upr <- ifelse(clim.bio$upr>1,1,clim.bio$upr)

save(clim.bio, file="clim_bio_all_variables.RData")


##### 2) Relative importance of climate and competition - from individual data #################################

load("clim_effect.RData")

ind.clim.bio.rel.imp <- as.data.frame(clim.effect[,c("latitude")])
names(ind.clim.bio.rel.imp) <- c("latitude")
ind.clim.bio.rel.imp$climate.rel <- NA
ind.clim.bio.rel.imp$competition.rel <- NA

# mean effect of the climate variables
clim.effect$effect_climate <- (abs(clim.effect$effect_PET) + abs(clim.effect$effect_minSPEI) +
                                       abs(clim.effect$effect_FROSTS)) /3


clim.effect$max <- apply(cbind(clim.effect$effect_climate, abs(clim.effect$effect_BAtotale)), 1, max)

ind.clim.bio.rel.imp$climate.rel <- clim.effect$effect_climate/clim.effect$max 
ind.clim.bio.rel.imp$competition.rel <- abs(clim.effect$effect_BAtotale)/clim.effect$max 

names(ind.clim.bio.rel.imp) <- c("Latitude", "Climate", "Competition")

save(ind.clim.bio.rel.imp, file="ind_rel_imp__climate_vs_compet.RData")
load("ind_rel_imp__climate_vs_compet.RData")

ind.clim.bio.rel.imp$Latitude_c <- cut(ind.clim.bio.rel.imp$Latitude, seq(36.5,70,0.5))

t.mean <- aggregate(ind.clim.bio.rel.imp[,c("Climate","Competition")], by=list(ind.clim.bio.rel.imp$Latitude_c), mean)
t.max <- aggregate(ind.clim.bio.rel.imp[,c("Climate","Competition")], by=list(ind.clim.bio.rel.imp$Latitude_c), max)
t.sd <- aggregate(ind.clim.bio.rel.imp[,c("Climate","Competition")], by=list(ind.clim.bio.rel.imp$Latitude_c), sd)
t.n <- aggregate(ind.clim.bio.rel.imp[,c("Climate","Competition")], by=list(ind.clim.bio.rel.imp$Latitude_c), length)

t.sd$climate.se <- t.sd$Climate/sqrt(t.n$Climate)
t.sd$competition.se <- t.sd$Competition/sqrt(t.n$Competition)

lats <- seq(36.5,69.5,0.5)
lats <- lats[!lats %in% c(37.5,38,38.5,39)]
t.mean$Latitude <- lats
t.sd$Latitude <- lats

ind.rel.imp.clim.bio.long <- melt(t.mean, id = "Latitude", measure = c("Climate","Competition"))
names(ind.rel.imp.clim.bio.long) <- c("Latitude", "variable", "mean")

## error bands
t.se <- t.sd[,4:6]
names(t.se) <- c("Climate","Competition", "Latitude")
ind.rel.imp.clim.bio.long.se <- melt(t.se, id = "Latitude", measure = c("Climate","Competition"))
names(ind.rel.imp.clim.bio.long.se) <- c("Latitude", "variable", "se")

# merge the mean and the se dataframes
clim.bio <- merge(ind.rel.imp.clim.bio.long, ind.rel.imp.clim.bio.long.se, by=c("Latitude", "variable"))
clim.bio$lwr <- clim.bio$mean-(1.96*clim.bio$se)
clim.bio$upr <- clim.bio$mean+(1.96*clim.bio$se) 

clim.bio$lwr <- ifelse(clim.bio$lwr<0,0,clim.bio$lwr)
clim.bio$upr <- ifelse(clim.bio$upr>1,1,clim.bio$upr)

save(clim.bio, file="clim_bio_climate_vs_compet.RData")
load("clim_bio_climate_vs_compet.RData")
