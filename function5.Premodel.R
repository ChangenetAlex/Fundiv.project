# Alex on the 06/06/2018
# Script to visualize my data before modeling and saving (or not) all the plots in one file 
# Function extraction + function Bootstrap

SavingInfo = "Premodel(
z = dfplot2, which is my dataframe with all variables I want to analyse
Resp=Resp, Which is my response variable, by defaut this is sp.mortality.plot.rate.yr
Explain=Explain, Which are the variable I want to study, by defaut a few are included 
size=4, For the cleveland plots, this is the par(mfrow=c(size,size)) parameter
save = T, to save in the new directory premodel or just display my plots)

Require the tidyr, lattice and usdm packages.

!!! WARNING !!! 
Before running this function you need to be in the directory corresponding to the SPECIES you want to study. 
!!! WARNING !!!

Example : Premodel(dfplot2,Resp=sp.mortality.ba,Explain=c(x,y...),size=3,save=T)

"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

## Packages ##

#library(usdm)
library(lattice)
library(tidyr)
library(usdm)

################### MY PARAMETERS ####################
Resp = c("sp.mortality.plot.count.yr")              ##
Explain = c("ppet.mean_climate_min.30",             ##
            "tmean.son_climate_max.30",             ##
            "BA.ha.plot.1",                         ##
            "BA.O.plot.1",                          ##
            "BAj.plot.1",                           ##
            "Plotcat")                              ##
######################################################

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

# Dotchart
Premodel <- function(z,Resp,Explain,size,save=F){
  Name <- 
        if(all(list.dirs(Dir,full.names=F)!="Premodel")) dir.create(paste0(Dir,"/Premodel/"))
        A <- cbind(Resp,Explain[-c(length(Explain))])
  c1 <- rep(z[,c(Resp)],times=length(Explain))
  c2 <- gather(z[,c(Explain)])
  test <- cbind(c1,c2[,c(2,1)])
  if (save==T) jpeg(file=paste0(Dir,"/Premodel/",Resp,"_Scatterplots.jpeg"),width=710)
  par(mfrow=c(1,1))
  xyplot(test[,1]~test[,2] | test[,3],col=1,
         strip=function(bg="white", ...)
           strip.default(bg="white", ...),
         scales=list(alternating=TRUE,
                     x=list(relation="free"),
                     y=list(relation="same")),
         xlab="explanatory variables",
         ylab="mortality",
         panel=function(x,y){
           panel.grid(h=-1,v=2)
           panel.points(x,y,col=1)
           panel.loess(x,y,xol=1,lwd=2)})
  dev.off()
  if (save==T) jpeg(file=paste0(Dir,"/Premodel/",Resp,"_ClevelandPlot.jpeg"),width=710) # Save it 
  par(mfrow=c(size,size))
  for (i in 1:length(A)){
    dotchart(z[,c(A[i])],
    groups=factor(z[,c(last(Explain))]),
    ylab = last(Explain),
    xlab = A[i],
    main = "Cleveland dotplot",
    pch = c(1:3))
  }
  dev.off()
  if (save==T) jpeg(file=paste0(Dir,"/Premodel/",Resp,"_",deparse(substitute(z)),"_",Explain[1],"_Hist.jpeg"),width=710) # Save it 
  par(mfrow=c(2,2))
  for (i in 2:length(A)){
    hist(z[,c(A[i])],breaks=40,xlab="Centered & normed Scale",ylab="Frequency",main="")
  }
  dev.off()
  if (save==T) jpeg(file=paste0(Dir,"/Premodel/",Resp,"_",deparse(substitute(z)),"_",Explain[1],"_PairsPlot.jpeg"),width=710) # Save it 
  par(mfrow=c(1,1))
  pairs(z[,c(A)],panel = panel.smooth,upper.panel=panel.cor,diag.panel = panel.hist)
  dev.off()
  VIF <- vif(z[,Explain[-c(length(Explain))]])
  VIFSTEP <- vifstep(z[,Explain[-c(length(Explain))]],th=8)
  VIF
  VIFSTEP
  if (save==T){
          capture.output(print(VIF),file=paste0(Dir,"/Premodel/",Resp,"_Vif.txt")) # Output as a latex wrapped in a txt file
          capture.output(print(VIFSTEP),file=paste0(Dir,"/Premodel/",Resp,"_Vifstep.txt")) # Output as a latex wrapped in a txt file
  }
}
