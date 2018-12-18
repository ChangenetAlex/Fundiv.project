library(knitr)
library(spaMM)
library(fields)
library(piecewiseSEM)
library(raster)
library(rworldmap)
library(rgdal)
library(rasterVis)
library(rgeos)
library(ggplot2)

i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allmod <- c("Mbin13A.18","Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin3B.20","Mbin15B.21","Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")
for (i in 1:length(Allcode)){
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  #rm("x")
  # New section 
  r1 <- data.frame(x$data$longitude, x$data$latitude, predict(x))
  colnames(r1) <- c("X","Y","Z")
  ?quantile
  if (x$family$family=="binomial") {r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"0",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"2","1"))} ## Percentage of mortality as classes 
  #if (x$family$family=="binomial") {r1$bin <- ifelse(r1$Z>0.20,"0",ifelse(r1$Z>0.1,"1","2"))} ## Percentage of mortality as classes 
  #}else r1$bin <- ifelse(r1$Z>40,"0",ifelse(r1$Z>20,"1","2"))
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(worldmap$REGION=="Europe"),]
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  europe = crop(europe,c(-15, 45, 35, 70))
  p <- ggplot() + theme(panel.background = element_rect(fill="lightblue", colour="black", size=3, 
                                                        linetype=1, color="black"),legend.key.size = unit(1.5, "cm"),legend.position = c(0.15,0.8)) +
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="black",size=0.5, linetype="solid", 
                                           colour ="black"))+
    geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill="gray20",col="gray14") + 
    coord_fixed(1.3) + geom_point(data = r1, aes(x = X, y = Y, group=bin, col = bin, shape = bin, size=bin))
  if (x$family$family=="binomial") {p <- p + scale_colour_manual(values = c("red","white","blue"),name="Probability of observing \n at least one event",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2)))) +
    scale_shape_manual(values = c(20,3,17),name="Probability of observing \n at least one event",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2)))) +
    scale_size_manual(values= c(0.7,0.5,0.7),name="Probability of observing \n at least one event",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2))))
  }else {p <- p + scale_colour_manual(values = c("red","white","blue"),name="Number of events by plot",labels = c(">40", "20-40","0-20")) +
    scale_shape_manual(values = c(20,3,17),name="Number of events by plot",labels = c(">40", "20-40","0-20")) +
    scale_size_manual(values= c(1,1,1),name="Number of events by plot",labels = c(">40", "20-40","0-20"))}
  p <- p + guides(shape = guide_legend(override.aes = list(size = 5)))+
    labs(title=paste0('Mortality occurence predicted probability of ',Allcode[i]), y=paste0("Latitude"), x="Longitude", caption="Changenet et al. 2018")+
    theme(text = element_text(face="bold"),legend.direction ="vertical",
          axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
          legend.key = element_rect(fill = "gray20", colour = "white"),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.5),
          plot.caption = element_text(face="bold.italic"))
  ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Figures.all/Synthesis/",Allcode[i],"_",Allmod[i],"_","Pred.Mort.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
}
  


i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allmod <- c("MnbZT13A.19","MnbZT14A.20","MnbZT13B.23","MnbZT13A.21","M2nbZT13B.29","M2nbZT13A.22","M2nbZT15B.24","MnbZT3C.24","M2nbZT7B.26","MnbZT11B.22","MnbZT3B.27","MnbZT15B.28","MnbZT5B.21","MnbZT13A.33","M2nbZT1C.27","M2nbZT1B.24","MnbZT13B.27","MnbZT7A.23","MnbZT1B.27")

for (i in 1:length(Allcode)){
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/Negbin/",Allmod[i],"/"))
  setwd(Dir)
  assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  #rm("x")
  # New section 
  r1 <- data.frame(x$data$longitude, x$data$latitude, predict(x))
  colnames(r1) <- c("X","Y","Z")
  r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"0",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"2","1")) ## Percentage of mortality as classes 
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(worldmap$REGION=="Europe"),]
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  europe = crop(europe,c(-15, 45, 35, 70))
  p <- ggplot() + theme(panel.background = element_rect(fill="lightblue", colour="black", size=3, 
                                                        linetype=1, color="black"),legend.key.size = unit(1.5, "cm"),legend.position = c(0.15,0.8)) +
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="black",size=0.5, linetype="solid", 
                                           colour ="black"))+
    geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill="gray20",col="gray14") + 
    coord_fixed(1.3) + geom_point(data = r1, aes(x = X, y = Y, group=bin, col = bin, shape = bin, size=bin))
  p <- p + scale_colour_manual(values = c("red","white","blue"),name="Rate of mortality events \n (proportion/plot/year)",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2)))) +
    scale_shape_manual(values = c(20,3,17),name="Rate of mortality events \n (proportion/plot/year)",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2)))) +
    scale_size_manual(values= c(0.9,0.5,0.9),name="Rate of mortality events \n (proportion/plot/year)",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2))))
  p <- p + guides(shape = guide_legend(override.aes = list(size = 5)))+
    labs(title=paste0('Mortality abundance predicted probability of ',Allcode[i]), y=paste0("Latitude"), x="Longitude", caption="Changenet et al. 2018")+
    theme(text = element_text(face="bold"),legend.direction ="vertical",
          axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
          legend.key = element_rect(fill = "gray20", colour = "white"),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.5),
          plot.caption = element_text(face="bold.italic"))
  ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Figures.all/Synthesis/",Allcode[i],"_",Allmod[i],"_","Pred.Mort.Abun.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
}



# Obs number

i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allmod <- c("MnbZT13A.19","MnbZT14A.20","MnbZT13B.23","MnbZT13A.21","M2nbZT13B.29","M2nbZT13A.22","M2nbZT15B.24","MnbZT3C.24","M2nbZT7B.26","MnbZT11B.22","MnbZT3B.27","MnbZT15B.28","MnbZT5B.21","MnbZT13A.33","M2nbZT1C.27","M2nbZT1B.24","MnbZT13B.27","MnbZT7A.23","MnbZT1B.27")
for (i in 1:length(Allcode)){
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/Negbin/",Allmod[i],"/"))
  setwd(Dir)
  assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  #rm("x")
  # New section 
  r1 <- data.frame(x$data$longitude, x$data$latitude, x$data$sp.mortality.plot.rate.yr)
  colnames(r1) <- c("X","Y","Z")
  r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"0",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"2","1")) ## Percentage of mortality as classes 
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(worldmap$REGION=="Europe"),]
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  europe = crop(europe,c(-15, 45, 35, 70))
  p <- ggplot() + theme(panel.background = element_rect(fill="lightblue", colour="black", size=3, 
                                                        linetype=1, color="black"),legend.key.size = unit(1.5, "cm"),legend.position = c(0.15,0.8)) +
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="black",size=0.5, linetype="solid", 
                                           colour ="black"))+
    geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill="gray20",col="gray14") + 
    coord_fixed(1.3) + geom_point(data = r1, aes(x = X, y = Y, group=bin, col = bin, shape = bin, size=bin))
  p <- p + scale_colour_manual(values = c("red","white","blue"),name="Proportion of observed \n mortality events/year",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2)))) +
    scale_shape_manual(values = c(20,3,17),name="Proportion of observed \n mortality events/year",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2)))) +
    scale_size_manual(values= c(0.9,0.5,0.9),name="Proportion of observed \n mortality events/year",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2))))
  p <- p + guides(shape = guide_legend(override.aes = list(size = 5)))+
    labs(title=paste0('Observed mortality of ',Allcode[i]), y=paste0("Latitude"), x="Longitude", caption="Changenet et al. 2018")+
    theme(text = element_text(face="bold"),legend.direction ="vertical",
          axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
          legend.key = element_rect(fill = "gray20", colour = "white"),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.5),
          plot.caption = element_text(face="bold.italic"))
  ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Figures.all/Synthesis/",Allcode[i],"_Mort.obs.Abun.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
}

# Obs number as a rate 


i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allmod <- c("Mbin13A.18","Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin3B.20","Mbin15B.21","Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")
for (i in 1:length(Allcode)){
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  #rm("x")
  # New section 
  r1 <- data.frame(x$data$longitude, x$data$latitude, x$data$sp.mortality.plot.rate.yr)
  colnames(r1) <- c("X","Y","Z")
  ?quantile
  if (x$family$family=="binomial") {r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"0",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"2","1"))} ## Percentage of mortality as classes 
  #if (x$family$family=="binomial") {r1$bin <- ifelse(r1$Z>0.20,"0",ifelse(r1$Z>0.1,"1","2"))} ## Percentage of mortality as classes 
  #}else r1$bin <- ifelse(r1$Z>40,"0",ifelse(r1$Z>20,"1","2"))
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(worldmap$REGION=="Europe"),]
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  europe = crop(europe,c(-15, 45, 35, 70))
  p <- ggplot() + theme(panel.background = element_rect(fill="lightblue", colour="black", size=3, 
                                                        linetype=1, color="black"),legend.key.size = unit(1.5, "cm"),legend.position = c(0.15,0.8)) +
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="black",size=0.5, linetype="solid", 
                                           colour ="black"))+
    geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill="gray20",col="gray14") + 
    coord_fixed(1.3) + geom_point(data = r1, aes(x = X, y = Y, group=bin, col = bin, shape = bin, size=bin))
  if (x$family$family=="binomial") {p <- p + scale_colour_manual(values = c("red","white","blue"),name="Probability of observing \n at least one event",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2)))) +
    scale_shape_manual(values = c(20,3,17),name="Probability of observing \n at least one event",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2)))) +
    scale_size_manual(values= c(0.7,0.5,0.7),name="Probability of observing \n at least one event",labels = c(paste0("> ",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0(signif(quantile(r1$Z,0.25,type=7),digits=2),"-",signif(quantile(r1$Z,0.75,type=7),digits=2)),paste0("< ",signif(quantile(r1$Z,0.25,type=7),digits=2))))
  }else {p <- p + scale_colour_manual(values = c("red","white","blue"),name="Number of events by plot",labels = c(">40", "20-40","0-20")) +
    scale_shape_manual(values = c(20,3,17),name="Number of events by plot",labels = c(">40", "20-40","0-20")) +
    scale_size_manual(values= c(1,1,1),name="Number of events by plot",labels = c(">40", "20-40","0-20"))}
  p <- p + guides(shape = guide_legend(override.aes = list(size = 5)))+
    labs(title=paste0('Mortality occurence predicted probability of ',Allcode[i]), y=paste0("Latitude"), x="Longitude", caption="Changenet et al. 2018")+
    theme(text = element_text(face="bold"),legend.direction ="vertical",
          axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
          legend.key = element_rect(fill = "gray20", colour = "white"),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.5),
          plot.caption = element_text(face="bold.italic"))
  ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Figures.all/Synthesis/",Allcode[i],"_",Allmod[i],"_","Pred.Mort.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
}



## Maps with no mortality plots and mortality observed as a continous variable 

i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allmod <- c("Mbin13A.18","Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin3B.20","Mbin15B.21","Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")
for (i in 1:length(Allcode)){
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  #rm("x")
  # New section 
  r1 <- data.frame(x$data$longitude, x$data$latitude, x$data$sp.mortality.plot.count.yr)
  r2 <- r1[r1$x.data.sp.mortality.plot.count.yr==0,]
  r1 <- r1[r1$x.data.sp.mortality.plot.count.yr>0,]
  r1$x.data.sp.mortality.plot.count.yr <- log(r1$x.data.sp.mortality.plot.count.yr)
  colnames(r1) <- c("X","Y","Z")
  colnames(r2) <- c("X","Y","Z")
  #r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"0",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"2","1")) ## Percentage of mortality as classes 
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(worldmap$REGION=="Europe"),]
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  europe = crop(europe,c(-15, 45, 35, 70))
  p <- ggplot() + theme(panel.background = element_rect(fill="white", size=3, 
                                                        linetype=1, color="black"),legend.key.size = unit(1.5, "cm"),legend.position = c(0.9,0.5)) +
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="white",size=0.1, linetype="solid", 
                                           colour ="black"))+
    geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill="gray70",col="gray70") + 
    coord_fixed(1.3) + 
    geom_point(data = r2, aes(x = X, y = Y),col="black",size=0.01,alpha=1) +
    #scale_discrete_manual(aes="fill") +
    geom_point(data = r1, aes(x = X, y = Y,col=Z,size=Z)) +
    scale_size(range = c(0.5,1.5),name="Mortality rate \n Events/year",guide=F) +
    scale_colour_gradient2(low="powderblue", mid="orange",high="firebrick4",midpoint=median(r1$Z),name="Mortality rate\n Events/year\n (log)")
  
  p
  p <-p + labs(title=paste0('Observed mortality (tree by years) of ',Allcode[i]), y=paste0("Latitude"), x="Longitude", caption="Changenet et al. 2018")+
    theme(text = element_text(face="bold"),legend.direction ="vertical",
          axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
          legend.key = element_rect(fill = "gray80", colour = "white"),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.5),
          plot.caption = element_text(face="bold.italic"))
  p
  ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Figures.all/Synthesis/",Allcode[i],"_Obs.Mort.Alltree.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
}

