# Alex on the 26/06/2018
# Adapted from Sophie Radcliffe and Juliette Archambau 
# Multiple functions to process the coeffficient of my models according to the latitude 

SavingInfo = "
####### function 1 ########

ExatracTest gives the parameters list of the model
Exemple : ExtracTest(Mbin1FAGSYLspaMM) # Me donne la liste des paramètres de mon modèles. 

####### function 2 ########

Effect_coef extract the calculation of the parameters that we want automatically (taking into account the log and the quadratic effect)
Example : 
for (i in c('BAIj.plot.bis','BA.ha.plot.1','BAj.plot.1','treeNbr','bio1_climate_mean.30','bio14_climate_mean.30')){
  Effect_coef(Mbin1FAGSYLspaMM,i)}

! ! !  WARNINGS  ! ! !  ! ! !  WARNINGS  ! ! !   ! ! !  WARNINGS  ! ! ! 

# In the model, interactions parameters have to be after log, quadratic and simple effects
# Also, transformed names have to be the same as non transformed names : For instance 
bio1_climate_mean.30 become logbio1_climate_mean.30 if logscaled. 

! ! !  WARNINGS  ! ! !  ! ! !  WARNINGS  ! ! !   ! ! !  WARNINGS  ! ! ! 


####### function 3 ########
Effect_summary is the Summary by latitude of the effects I want and effects regroupement (competition,biotic,climate) # Cf script for more details 

Example : Effect_summary(x,BioticPara='competition')

####### function 4 ########
ggeffect allow to plot this parameters
x is my model
y is the kind of info that i want either 'ABS' for absolute values (compared to 0), either 'REL' for comparison between each parameters
effects is the parameters we xould like to plot. it can be either 'indiv' for simple effects, or 'sum' for synthetic effect such as copetition or climate

Example : ggEffect <- function(x,y='REL',effect='sum')
"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))



################################################################################

################################################################################


#### Load libraries #######
library(reshape2)         #
require(raster)           #
library(ggplot2)          #
library(mgcv)             # 
library(grid)             #
#library(Cairo)           #
###########################



### Obtain the para #######

ExtracTest <- function(x){
  A <- paste0(getCall(x)[2])
  A <- unlist(strsplit(A, "~",fixed = T))[2]
  A <- unlist(strsplit(A, "+",fixed = T))
  A <- grep(A,pattern = "|",fixed=T,value=T,invert=T)
  A <- grep(A,pattern = ":",fixed=T,value=T,invert=T)
  A <- grep(A,pattern = "e)",fixed=T,value=T,invert=T) #Remove the AC term which remains
  A <- grep(A,pattern = "^[ ]$",value=T,invert=T) #Remove empty characters
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


### Calculation of each effect that I specify ###

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

### Summary by latitude of the effects I want and effects regroupement (competition,biotic, climate) ###

Effect_summary <- function(x,BioticPara="competition"){
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
  if (BioticPara != "competition") {Var.biotic <- grep(Var,pattern = paste(c("BA","tree","dbh"),collapse="|"),value=T,invert=F)     #                                    #
  #################                                                                                                           #
  ## Competition ##    # Same as Biotic But the BAI is removed since it is a growth param                                     #
  #################                                                                                                           #
  } else {Var.comp <- grep(Var,pattern = paste(c("BA","tree"),collapse="|"),value=T,invert=F)                                 #
  Var.comp <- grep(Var.comp,pattern = "BAI",fixed=T,value=T,invert=T)}                                                      #
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
    ind.clim.bio.rel.imp[,EffectCol[i]] <- abs(clim.effect[,EffectCol[i]])/clim.effect$max  # Importance compare to the max one tha was extracted
    #ind.clim.bio.rel.imp[,EffectCol[i]] <- abs(clim.effect[,EffectCol[i]])/clim.effect$sum  # Relative importance compare to the others (sum is qual to one)
    ind.clim.bio.abs.imp[,EffectCol[i]] <- abs(clim.effect[,EffectCol[i]])  # Relative importance compare to the others (sum is qual to one)
  }
  for (i in 1:length(EffectCum)){
    ind.clim.bio.rel.imp[,EffectCum[i]] <- abs(clim.effect[,EffectCum[i]])/clim.effect$maxcum  # Importance compare to the max one
    #ind.clim.bio.rel.imp[,EffectCum[i]] <- abs(clim.effect[,EffectCum[i]])/clim.effect$sumcum  # Relative importance compare to the others (sum is qual to one)
    ind.clim.bio.abs.imp[,EffectCum[i]] <- abs(clim.effect[,EffectCum[i]])  # Relative importance compare to the others (sum is qual to one)
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

### Representation on a plot of the parameters extracted 
ggEffect <- function(x,y="REL",effect="sum",band=T){
  if (y=="ABS"){clim.bio <- get(load(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",
                                            CODE,"/CLIMAP/Models/binomial/",
                                            deparse(substitute(x)),"/clim_bio_abs_",
                                            deparse(substitute(x)),".RData")))
  }else if(y=="REL"){clim.bio <- get(load(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",
                                                 CODE,"/CLIMAP/Models/binomial/",
                                                 deparse(substitute(x)),"/clim_bio_rel_",
                                                 deparse(substitute(x)),".RData")))}
  ## Effect names 
  EffectCol <- grep(as.character(unique(clim.bio$variable)),pattern = "effect",fixed=T,value=T,invert=F) # Simple effect columns for which coef were extracted
  
  ## Individual effects vs sum of the effects 
  if (effect=="indiv"){clim.bio <- clim.bio[clim.bio$variable%in%EffectCol,] # If individual, keep these individual effects 
  }else if (effect=="sum"){
    clim.bio <- clim.bio[!(clim.bio$variable%in%EffectCol),] # If not, keep all the orther effects who have not "effect in their names"
    EffectCol <- as.character(unique(clim.bio$variable))}    # Exatrct these sum effects and put it in the initial vector
  clim.bio[,"variable"] <- as.character(clim.bio[,"variable"]) # Transfrom the factor as a charcater vector
  
  # Find my missing values 
  A <- unique(clim.bio[is.na(clim.bio$se),"latitude"])
  missing <- data.frame(xmin=A,xmax=A+0.5, ymin=rep(Inf,length=length(A)), ymax=rep(-Inf,length=length(A))) # Identify the missing values automatically 
  i = 1
  while (is.na(missing$xmin[i+1])==F){
    if (missing$xmax[i+1]-missing$xmax[i]==0.5){
      while (missing$xmax[i+1]-missing$xmax[i]==0.5 & is.na(missing$xmax[i+1])==F){
        missing$xmax[i] <- missing$xmax[i+1]
        missing <- missing[-c(i+1),]
      }} else i=i+1
  }
  
  Mycol <- c("red","dark green", "dodgerblue3","orange","yellow","gray","black","green","lightblue","gray20") #### The colors that I will use 
  Myline <- c(1,3,6,1,3,6,1,3,6,1,3,6)
  
  # The plot theme 
  p<-ggplot(clim.bio) + 
    theme_bw()
  if (y=="ABS"){p <- p + ylab("Predicted absolute importance") +
    scale_y_continuous(limits=c(-0.2,max(clim.bio$mean)))
  }else if(y=='REL'){p <- p + ylab("Predicted relative importance") + 
    scale_y_continuous(limits=c(-0.2,1))}
  p <- p + theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
  ) +
    theme(axis.text.x = element_text(size=11),
          text = element_text(face="bold"),#
          legend.background=element_rect(fill="white",colour="black",size=0.2),#
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),#
          axis.line = element_line(colour="black"),#
          plot.title = element_text(size=14,hjust = 0.5),#
          plot.caption = element_text(face="bold.italic"),#
          axis.text.y = element_text(size=11),  
          axis.title.x = element_text(size=11),
          axis.title.y = element_text(size=11),
          legend.position="bottom",
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.key.size = unit(5, 'lines'),
          legend.text = element_text(size=11,face="bold"),
          legend.key.height=unit(2,"line"),
          legend.key.width=unit(5,"line")
    )
  # The plot 
  par(mar=c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  p <- p + geom_line(aes(latitude, mean, colour = as.character(variable), linetype=as.character(variable)), size=0.8)+
    guides(fill=FALSE) +
    guides(col=guide_legend(ncol=3, byrow=F))
  if (effect=="indiv"){p <- p + scale_colour_manual(values=c(Mycol[1:length(EffectCol)]),labels=paste0("  ",EffectCol,"  ")) +
    scale_linetype_manual(values=c(Myline[1:length(EffectCol)]),labels=paste0("  ",EffectCol,"  "))
  } else if (effect=="sum"){p <- p + scale_colour_manual(values=c("darkorchid4","orange"),labels=paste0("  ",EffectCol,"  "))}
  
  if (band==T){p <- p + geom_ribbon(data=clim.bio,aes(latitude, mean, ymin=lwr, ymax=upr, colour=variable, fill=variable),alpha=0.05, linetype=2)
  }else p <- p
  p <- p + geom_rect(data=missing, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),colour="light grey", fill="white", inherit.aes=FALSE) +
    geom_rect(xmin=min(clim.bio$latitude), xmax=max(clim.bio$latitude), ymin=-0.10, ymax=-0.05,colour="black", fill="black")+
    geom_rect(xmin=0,xmax=42.5,ymin=-Inf,ymax=-0.1,colour="black", fill="red",alpha=0.2)+
    geom_rect(xmin=42.5,xmax=58,ymin=-Inf,ymax=-0.1,colour="black", fill="green",alpha=0.5)+
    geom_rect(xmin=58,xmax=Inf,ymin=-Inf,ymax=-0.1,colour="black", fill="blue",alpha=0.8)+
    geom_label(label = "Mediterranean", x=42.5, y=-0.01, size=3.5,vjust=2.2,label.padding = unit(0.15, "lines"),
               label.r = unit(0.15, "lines"),colour="white",fill="red",label.size = 0.2)+
    geom_label(label = "Temperate", x=50, y=-0.01, size=3.5,vjust=2.2,label.padding = unit(0.15, "lines"),
               label.r = unit(0.15, "lines"),colour="white",fill="darkgreen",label.size = 0.2)+
    geom_label(label = "Boreal", x=58, y=-0.01, size=3.5,vjust=2.2,label.padding = unit(0.15, "lines"),
               label.r = unit(0.15, "lines"),colour="white",fill="blue",label.size = 0.2)+
    geom_segment(x=min(clim.bio$latitude), y=-0.1, xend=min(clim.bio$latitude), yend=Inf, colour="black", size=0.1,linetype=11) +
    geom_segment(x=max(clim.bio$latitude), y=-0.1, xend=max(clim.bio$latitude), yend=Inf, colour="black", size=0.1,linetype=11) +
    geom_segment(x=58, y=-0.1, xend=58, yend=Inf, colour="light grey", size=0.1,linetype=2) +
    geom_segment(x=42.5, y=-0.1, xend=42.5,yend=Inf, colour="light grey", size=0.1,linetype=2) +
    labs(caption="Changenet et al. 2018")
  print(p)
  if (band==T){ggsave(filename = paste0(y,"_",effect,"_",CODE,"_",deparse(substitute(x)),"_band.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in") # Save 
  }else ggsave(filename = paste0(y,"_",effect,"_",CODE,"_",deparse(substitute(x)),"_NO.band.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
}

# width = 12.37, height = 7.04, dpi=150,units = "in"


