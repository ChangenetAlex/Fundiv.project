Mbin1FAGSYLspaMM <- fitme(sp.mort.bin ~ 0 + treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1 + bio14_climate_mean.30 + logbio14 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1 + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:logbio14 + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:logbio14 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1:logbio14 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:logbio1 + Plotcat:logbio14 + (1|country), data=dfplot2,family = binomial(),method='REML')
Mbin1FAGSYLspaMM <- fitme(sp.mort.bin ~ treeNbr + yearsbetweensurveys + min_spei01 + mean_spei01 + min_spei06 + mean_spei06 + min_spei12 + mean_spei12 + min_spei24 + mean_spei24 + min_spei48 + mean_spei48 + BAIj.plot.bis + BA.ha.plot.1 + BAj.plot.1 + bio1_climate_mean.30 + logbio1_climate_mean.30 + bio14_climate_mean.30 + logbio14_climate_mean.30 + Plotcat + I(bio1_climate_mean.30^2) + I(bio14_climate_mean.30^2) + BA.ha.plot.1:bio1_climate_mean.30 + BA.ha.plot.1:logbio1_climate_mean.30 + BAj.plot.1:bio1_climate_mean.30 + BAj.plot.1:logbio1_climate_mean.30 + BA.ha.plot.1:bio14_climate_mean.30 + BA.ha.plot.1:logbio14_climate_mean.30 + BAj.plot.1:bio14_climate_mean.30 + BAj.plot.1:logbio14_climate_mean.30 + BA.ha.plot.1:BAj.plot.1 + bio1_climate_mean.30:bio14_climate_mean.30 + logbio1_climate_mean.30:logbio14_climate_mean.30 + BA.ha.plot.1:Plotcat + BAj.plot.1:Plotcat + Plotcat:bio1_climate_mean.30 + Plotcat:bio14_climate_mean.30 + Plotcat:logbio1_climate_mean.30 + Plotcat:logbio14_climate_mean.30 + (1|country), data=dfplot2,family = binomial(),method='REML')
Saving(Mbin1FAGSYLspaMM)


# Dans mon modèle, les paramètres d'interactions doivent être après les para log, seuls et cubique
# Les noms des variables scaled doivent être les même que les non scales et pas une abbréviation

# get the parameter values from the glm model Unique parameters of the model
library(reshape2)
require(raster)
library(ggplot2)
library(mgcv)
library(grid)
library(Cairo)

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
Effect_summary <- function(x,BioticPara){
  if (grepl(deparse(substitute(x)),pattern="bin",fixed=T)==T){
    Dir <- paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/FAGSYL/CLIMAP/Models/binomial/",deparse(substitute(x)),"/")
  }else Dir <- paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/FAGSYL/CLIMAP/Models/Negbin/",deparse(substitute(x)),"/")
  setwd(Dir)      
  clim.effect <- get(load(paste0("clim_effect_",deparse(substitute(x)),".RData"))) ## Load my effects coefficients database (product from the second function)
  
  ############################################
  #### All effects that we will summarize ####
  ############################################
  
  EffectCol <- grep(colnames(clim.effect),pattern = "effect",fixed=T,value=T,invert=F) # Simple effect columns for which coef were extracted
  EffectCum <- c("climate","biotic","competition") # Add three combined categories
  if (BioticPara == "competition") EffectCum <- EffectCum[c(1,3)] else EffectCum <- EffectCum[1:2]
  EffectAll <- c(EffectCol,EffectCum) ## Simple + Categories 
  # Newdataframes made of these extracted effect and the latitudes 
  ind.clim.bio.abs.imp <- as.data.frame(clim.effect[,c("latitude")]) ## One with absolute values
  names(ind.clim.bio.abs.imp) <- c("latitude")
  ind.clim.bio.rel.imp <- as.data.frame(clim.effect[,c("latitude")]) ## The other with relative values
  names(ind.clim.bio.rel.imp) <- c("latitude")
  
  # Create as many columns as the number of simple effects calculted for both df, and fil it with NA
  for (i in 1:length(EffectCol)){
    ind.clim.bio.rel.imp[,EffectCol[i]] <- NA
    ind.clim.bio.abs.imp[,EffectCol[i]] <- NA
  }
  
  for (i in 1:length(EffectCum)){
    ind.clim.bio.rel.imp[,EffectCum[i]] <- NA
    ind.clim.bio.abs.imp[,EffectCum[i]] <- NA
  }
  
  #############################################################################################################################
  ######                                                                                                           ############
  ###### In the next section we define the paramters to include in the summary of biotic, competition and climatic ############
  ######                  To change them, we need to modify directly the 'paste0("...")' terms                     ############
  ######                                                                                                           ############
  #############################################################################################################################
  #
  Var <- grep(colnames(clim.effect),pattern = "effect",fixed=T,value=T,invert=F) # Effect terms, no log nor quadratique effect# 
  ############                                                                                                                #
  ## Biotic ##      ##  'BA' terms and the number of trees (if an effect was calculated for them)n + BAI                      #
  ############                                                                                                                #
  if (BioticPara != "competition") {Var.biotic <- grep(Var,pattern = paste(c("BA","tree"),collapse="|"),value=T,invert=F)     #                                    #
  #################                                                                                                           #
  ## Competition ##    # Same as Biotic But the BAI is removed since it is a growth param                                     #
  #################                                                                                                           #
  } else {Var.comp <- grep(Var,pattern = paste(c("BA","tree"),collapse="|"),value=T,invert=F)                                 #
  Var.comp <- grep(Var.biotic,pattern = "BAI",fixed=T,value=T,invert=T)}                                                      #
  ################                                                                                                            #
  ##  Climatic  ##     ## Variable de climat (spei + worlclim)                                                                #
  ################                                                                                                            #
  Var.climate <- grep(Var,pattern = paste(c("spei","climate"),collapse="|"),value=T,invert=F)                                 #
  #
  #############################################################################################################################
  ######                                                                                                           ############
  ######                                  This is the end of this section                                          ############
  ######                                                                                                           ############
  #############################################################################################################################
  
  
  # For each latitude, the sum of the effects corresponding to the three defined categories is calculated as a new effect column 
  if (BioticPara == "competition"){
    if (length(Var.comp)>1){clim.effect$competition <- apply(abs(clim.effect[,Var.comp]), 1, function(x) sum(x)/length(Var.comp)) #need length superior to 1
    }else clim.effect$competition <- abs(clim.effect[,Var.comp])
    if (length(Var.climate)>1){clim.effect$climate <- apply(abs(clim.effect[,Var.climate]), 1, function(x) sum(x)/length(Var.climate)) #need length superior to 1
    }else clim.effect$climate <- abs(clim.effect[,Var.climate])
  }else{
    if (length(Var.biotic)>1){clim.effect$biotic <- apply(abs(clim.effect[,Var.biotic]), 1, function(x) sum(x)/length(Var.biotic)) #need length superior to 1
    }else clim.effect$biotic <- abs(clim.effect[,Var.biotic])
    if (length(Var.climate)>1){clim.effect$climate <- apply(abs(clim.effect[,Var.climate]), 1, function(x) sum(x)/length(Var.climate)) #need length superior to 1
    }else clim.effect$climate <- abs(clim.effect[,Var.climate])}
  
  ## Exatrct abs and relative effects for single or categories
  clim.effect$max <- apply(cbind(abs(clim.effect[,EffectCol])), 1, max) # Extract Maximum among the simple effects
  clim.effect$sum <- apply(cbind(abs(clim.effect[,EffectCol])), 1, sum) # Sum of the simple effects
  clim.effect$maxcum <- apply(cbind(abs(clim.effect[,EffectCum])), 1, max) # Max between Competition and climate
  clim.effect$sumcum <- apply(cbind(abs(clim.effect[,EffectCum])), 1, sum) # Sum of competition and climate 
  
  for (i in 1:length(EffectCol)){
    ind.clim.bio.abs.imp[,EffectCol[i]] <- abs(clim.effect[,EffectCol[i]])/clim.effect$max  # Importance compare to the max one tha was extracted
    ind.clim.bio.rel.imp[,EffectCol[i]] <- abs(clim.effect[,EffectCol[i]])/clim.effect$sum  # Relative importance compare to the others (sum is qual to one)
  }
  for (i in 1:length(EffectCum)){
    ind.clim.bio.abs.imp[,EffectCum[i]] <- abs(clim.effect[,EffectCum[i]])/clim.effect$maxcum  # Importance compare to the max one
    ind.clim.bio.rel.imp[,EffectCum[i]] <- abs(clim.effect[,EffectCum[i]])/clim.effect$sumcum  # Relative importance compare to the others (sum is qual to one)
  }
  ## Save these two files 
  save(ind.clim.bio.abs.imp, file="./ind_abs_imp_all_variables.RData") # Save these two files
  save(ind.clim.bio.rel.imp, file="./ind_rel_imp_all_variables.RData")
  
  ### Create a synthetic information with cut latitudes
  
  ## Cut latitude in categories (min and max of my data) by sequences of 0.5°
  ind.clim.bio.rel.imp$latitude_c <- cut(ind.clim.bio.rel.imp$latitude, seq(floor(min(ind.clim.bio.rel.imp$latitude)),ceiling(max(ind.clim.bio.rel.imp$latitude)),0.5))
  ind.clim.bio.abs.imp$latitude_c <- cut(ind.clim.bio.abs.imp$latitude, seq(floor(min(ind.clim.bio.abs.imp$latitude)),ceiling(max(ind.clim.bio.abs.imp$latitude)),0.5))
  
  ## For simple and categories effects : calculate max, number, mean and standard deviation BY GROUP OF LATITUDE
  t.mean.rel <- ddply(ind.clim.bio.rel.imp, .(latitude_c), colwise(mean,c(EffectAll)), .drop=F,na.rm=T) #ok 
  #t.max.rel <- ddply(ind.clim.bio.rel.imp, .(latitude_c), colwise(max,na.rm=T,c(EffectAll)),.drop=F,na.rm=T,.inform = T) # Pbm with NA
  t.sd.rel <- ddply(ind.clim.bio.rel.imp, .(latitude_c), colwise(sd,c(EffectAll)), .drop=F,na.rm=T) #ok
  t.n.rel <- ddply(ind.clim.bio.rel.imp, .(latitude_c), colwise(length,c(EffectAll)), .drop=F) #ok
  
  t.mean.abs <- ddply(ind.clim.bio.abs.imp, .(latitude_c), colwise(mean,c(EffectAll)), .drop=F,na.rm=T) #ok 
  #t.max.abs <- ddply(ind.clim.bio.abs.imp, .(latitude_c), colwise(max,na.rm=T,c(EffectAll)),.drop=F,na.rm=T,.inform = T) # Pbm with NA
  t.sd.abs <- ddply(ind.clim.bio.abs.imp, .(latitude_c), colwise(sd,c(EffectAll)), .drop=F,na.rm=T) #ok
  t.n.abs <- ddply(ind.clim.bio.abs.imp, .(latitude_c), colwise(length,c(EffectAll)), .drop=F) #ok
  
  ## Calculation of the standard error from sd and number, by categories of latitude
  for (i in 1:length(EffectAll)){
    t.sd.rel[,paste0(EffectAll[i],"_se")] <- t.sd.rel[,EffectAll[i]]/sqrt(t.n.rel[,EffectAll[i]]) # Calcul of standard errors 
    t.sd.abs[,paste0(EffectAll[i],"_se")] <- t.sd.abs[,EffectAll[i]]/sqrt(t.n.abs[,EffectAll[i]]) # Calcul of standard errors 
  }
  
  ## label of latitude (0.5 superior degree)
  lats <- seq(floor(min(ind.clim.bio.rel.imp$latitude)),ceiling(max(ind.clim.bio.rel.imp$latitude)),0.5)[-1] # Remove the first 
  t.mean.rel$latitude <- lats
  t.sd.rel$latitude <- lats
  t.mean.abs$latitude <- lats
  t.sd.abs$latitude <- lats
  
  ## shape my mean data for ggplots (abs and relative)
  ind.rel.imp.clim.bio.long <- melt(t.mean.rel, id = "latitude", measure = c(EffectAll))
  names(ind.rel.imp.clim.bio.long) <- c("latitude", "variable", "mean")
  ind.abs.imp.clim.bio.long <- melt(t.mean.abs, id = "latitude", measure = c(EffectAll))
  names(ind.abs.imp.clim.bio.long) <- c("latitude", "variable", "mean")
  ## error bands for ggplots (se)
  t.se.rel <- t.sd.rel[,(ncol(t.sd.rel)-length(EffectAll)):ncol(t.sd.rel)]
  colnames(t.se.rel) <- sub("_se","", colnames(t.se.rel), ignore.case = FALSE,fixed = T) # remove the SE in the name to match all the variable
  ind.rel.imp.clim.bio.long.se <- melt(t.se.rel, id = "latitude", measure = c(EffectAll))
  names(ind.rel.imp.clim.bio.long.se) <- c("latitude", "variable", "se")
  t.se.abs <- t.sd.abs[,(ncol(t.sd.abs)-length(EffectAll)):ncol(t.sd.abs)] #same thing for relative 
  colnames(t.se.abs) <- sub("_se","", colnames(t.se.abs), ignore.case = FALSE,fixed = T)
  ind.abs.imp.clim.bio.long.se <- melt(t.se.abs, id = "latitude", measure = c(EffectAll))
  names(ind.abs.imp.clim.bio.long.se) <- c("latitude", "variable", "se")
  
  # Merge the mean and the se dataframes
  clim.bio.rel <- merge(ind.rel.imp.clim.bio.long, ind.rel.imp.clim.bio.long.se, by=c("latitude", "variable"))
  clim.bio.rel$lwr <- clim.bio.rel$mean-(1.96*clim.bio.rel$se) # lower intervalle
  clim.bio.rel$upr <- clim.bio.rel$mean+(1.96*clim.bio.rel$se) # higher intervalle 
  clim.bio.rel$lwr <- ifelse(clim.bio.rel$lwr<0,0,clim.bio.rel$lwr) 
  clim.bio.rel$upr <- ifelse(clim.bio.rel$upr>1,1,clim.bio.rel$upr)
  
  clim.bio.abs <- merge(ind.abs.imp.clim.bio.long, ind.abs.imp.clim.bio.long.se, by=c("latitude", "variable"))
  clim.bio.abs$lwr <- clim.bio.abs$mean-(1.96*clim.bio.abs$se)
  clim.bio.abs$upr <- clim.bio.abs$mean+(1.96*clim.bio.abs$se) 
  clim.bio.abs$lwr <- ifelse(clim.bio.abs$lwr<0,0,clim.bio.abs$lwr)
  clim.bio.abs$upr <- ifelse(clim.bio.abs$upr>1,1,clim.bio.abs$upr)
  
  ## Save it all 
  save(clim.bio.rel, file=paste0("clim_bio_rel_",deparse(substitute(x)),".RData"))
  save(clim.bio.abs, file=paste0("clim_bio_abs_",deparse(substitute(x)),".RData"))
}


ExtracTest(Mbin1FAGSYLspaMM) # Me donne la liste des paramètres de mon modèles. 
for (i in c("BAIj.plot.bis","BA.ha.plot.1","BAj.plot.1","treeNbr","bio1_climate_mean.30","bio14_climate_mean.30")){
  Effect_coef(Mbin1FAGSYLspaMM,i)}     # Para que l'on veut regarder   
Effect_summary(Mbin1FAGSYLspaMM,"competition")


### ggplot part = 4eme function  # to write 
# un para is x (pour la base de données) et un y pour le effect individuel VS effet sum 
# et un para pour les bandes ou pas 

clim.bio.abs <- clim.bio.abs[clim.bio.abs$variable%in%EffectCol,] # If individual effect, this is what we want ! 
clim.bio.abs <- clim.bio.abs[!(clim.bio.abs$variable%in%EffectCol),] # If only sum effect 
clim.bio.abs[,"variable"] <- as.character(clim.bio.abs[,"variable"])

p<-ggplot(clim.bio.abs) + 
  theme_bw() + 
  #theme_light(base_size = 15)+
  ylab("Predicted relative importance") + 
  scale_y_continuous(limits=c(-0.01,1)) +

  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
  ) +
  
  theme(axis.text.x = element_text(size=13),
        text = element_text(face="bold"),#
        legend.background=element_rect(fill="white",colour="black",size=0.2),#
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),#
        axis.line = element_line(colour="black"),#
        plot.title = element_text(size=18,hjust = 0.5),#
        plot.caption = element_text(face="bold.italic"),#
        axis.text.y = element_text(size=12),  
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(size=12),
        legend.key.height=unit(2,"line"),
        legend.key.width=unit(4,"line")
  )


missing <- data.frame(xmin=55,xmax=56.5, ymin=Inf, ymax=-Inf) # Identify the missing values automatically 
Mycol <- c("red","dark green", "dodgerblue3","orange","yellow","gray")
EffectCol <- unique(clim.bio.abs[,"variable"])
par(mar=c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
p <- p + geom_line(aes(latitude, mean, colour = as.character(variable)), size=0.8) +
  guides(linetype=FALSE, fill=FALSE) +
  guides(col=guide_legend(ncol=3, byrow=F)) +
  scale_colour_manual(values=c(Mycol[1:length(EffectCol)]),
                      labels=EffectCol) +
  geom_ribbon(data=clim.bio.abs,aes(latitude, mean, ymin=lwr, ymax=upr, colour=variable, fill=variable),alpha=0.05, linetype=2)+
   geom_rect(data=missing, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),colour="light grey", fill="white", inherit.aes=FALSE) +
  geom_rect(xmin=min(clim.bio.abs$latitude), xmax=max(clim.bio.abs$latitude), ymin=-0.025, ymax=0,colour="black", fill="black")+
  geom_rect(xmin=0,xmax=42.5,ymin=-Inf,ymax=-0.025,colour="black", fill="red",alpha=0.2)+
  #annotate("text", label = "Mediterranean", x=39.25, y=-0.01, vjust=1.2, size=3.3)+
  geom_rect(xmin=42.5,xmax=58,ymin=-Inf,ymax=-0.025,colour="black", fill="green",alpha=0.5)+
  geom_rect(xmin=58,xmax=Inf,ymin=-Inf,ymax=-0.025,colour="black", fill="blue",alpha=0.8)+
    geom_text(label = "Mediterranean", x=42.5, y=-0.01, size=4,vjust=2.5) +
    geom_text(label = "Temperate", x=50, y=-0.01, size=4,vjust=2.5) +
    geom_text(label = "Boreal", x=58, y=-0.01, size=4,vjust=2.5) +
  geom_segment(x=min(clim.bio.abs$latitude), y=-0.02, xend=min(clim.bio.abs$latitude), yend=Inf, colour="black", size=0.1,linetype=11) +
  geom_segment(x=max(clim.bio.abs$latitude), y=-0.02, xend=max(clim.bio.abs$latitude), yend=Inf, colour="black", size=0.1,linetype=11) +
  geom_segment(x=58, y=-0.02, xend=58, yend=Inf, colour="light grey", size=0.1,linetype=2) +
  geom_segment(x=42.5, y=-0.02, xend=42.5,yend=Inf, colour="light grey", size=0.1,linetype=2) +
  labs(caption="Changenet et al. 2018")
p


print(p)
dev.off()
