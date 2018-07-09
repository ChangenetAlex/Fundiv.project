# Alex on the 01/06/2018
# Script to compute the confidence interval with spaMM for any interaction 
# Function extraction + function Bootstrap

SavingInfo = " This function is to make the diagnostic and the cross validation for the mortality models 
Diagnostic(
x = Mymodel
Yportion = Percentage of the data on which the cross validation is done
AllInOne = Do you want all the figure to be on a single one or in different ones 
)
!!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! 

Before running this function you need to be in the directory corresponding to the model you want to study. 

!!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! 

Output:
Several plot and the a txt file with the values of the cross validation
"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

## Packages 

## Packages
library(spaMM)
library(dplyr)
library(arm)

## The function 
x <- M2bin15B
Diagnostic <- function(x,Yportion,AllInOne=T){
        
        # Extract the response and the explicative variable 
        ExplVar <- paste0(getCall(x)[2])
        Response <- unlist(strsplit(ExplVar, " ~ ",fixed = T))[1]
        ExplVar <- unlist(strsplit(ExplVar, " ~ ",fixed = T))[2]
        ExplVar <- unlist(strsplit(ExplVar, "+",fixed = T))
        ExplVar <- grep(ExplVar,pattern = "|",fixed=T,value=T,invert=T)
        ExplVar <- grep(ExplVar,pattern = ":",fixed=T,value=T,invert=T)
        ExplVar <- grep(ExplVar,pattern = "e)",fixed=T,value=T,invert=T) #Remove the AC term which remains
        ExplVar <- grep(ExplVar,pattern = "^[ ]$",value=T,invert=T) #Remove empty characters
        ExplVar <- sub("I(","", ExplVar, ignore.case = FALSE,fixed = T)
        ExplVar <- sub("^2)","", ExplVar, ignore.case = FALSE,fixed = T)
        ExplVar <- sub("offset(log(","", ExplVar, ignore.case = FALSE,fixed = T)
        ExplVar <- sub("))","", ExplVar, ignore.case = FALSE,fixed = T)
        ExplVar <- sub("\n ","", ExplVar, ignore.case = FALSE,fixed = T) # New line
        ExplVar <- sub(" ","", ExplVar, ignore.case = FALSE,fixed = T)
        ExplVar <- sub(" ","", ExplVar, ignore.case = FALSE,fixed = T) #Twice
        ExplVar <- unique(ExplVar)
        
        ### Cross validation ###
        
        if (x$family$family=="negbin"){
                setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/Negbin/",deparse(substitute(x)),"/"))
                #### Define my data ##### That are all data I have not taken in my model
                calibrate.data=sample(1:nrow(dfplot2[dfplot2$sp.mortality.plot.count>0,]), Yportion*nrow(dfplot2[dfplot2$sp.mortality.plot.count>0,])) # 66% of our data 
                df1=dfplot2[!is.na(dfplot2$BAI.O.plot.1)&dfplot2$sp.mortality.plot.count>0,][calibrate.data,]
                df2=dfplot2[!is.na(dfplot2$BAI.O.plot.1)&dfplot2$sp.mortality.plot.count>0,][-calibrate.data,]
          
          
                M1 = update(x, data=df1)
                M2 = predict (M1, newdata=df2,re.form=NA,type="response")
                
                ### Cor Test ###
                if (AllInOne == T){jpeg(file=paste0(deparse(substitute(x)),"Diagnostic.jpeg"),width=710);par(mfrow = c(2,2))
                } else {jpeg(file=paste0(deparse(substitute(x)),"_Response_VS_predicted.jpeg"),width=710)}
                D <- cor.test(df2$sp.mortality.plot.count.yr, M2, method = ("pearson")) # To record and store somewhere
                capture.output(print(D), file=paste0(deparse(substitute(x)),"_CrossValid.Table.txt")) # Output as a latex wrapped in a txt file
                plot(df2$sp.mortality.plot.count.yr, M2,xlab="Fitted values", ylab="Predicted", main="Predicted vs. fitted",
                     cex.main=1.5, cex.lab =1.3) # To save 
                abline(0,1)
                if (AllInOne == F){dev.off()}
                
                ### Residuals ###
                ResM <- residuals(M1)
                FitM <- fitted(M1)
                if (AllInOne == F) jpeg(file=paste0(deparse(substitute(x)),"_Residuals_VS_fitted.jpeg"),width=710) # Save it 
                plot(ResM ~ FitM, xlab="Fitted values", ylab="Residuals", main="Residuals vs. fitted",
                     cex.main=1.5, cex.lab =1.3)
                abline(h=0)
                if (AllInOne == F){dev.off()} # Save it 
                
                # Save this plot 
                if (AllInOne == F){jpeg(file=paste0(deparse(substitute(x)),"_QQplot.jpeg"),width=710)} # Save it 
                qqnorm(ResM, cex.main=1.5, cex.lab =1.3)
                qqline(ResM)
                if (AllInOne == F){dev.off()} # Save it 
                
                # Save this one 
                if (AllInOne == F){jpeg(file=paste0(deparse(substitute(x)),"_Observed_VS_fitted.jpeg"),width=710)} # Save it 
                plot(M1$data$sp.mortality.plot.count.yr ~ FitM, xlab="Fitted values", ylab="Observed", main="Observed vs. fitted",
                     cex.main=1.5, cex.lab =1.3)
                abline(0,1)
                dev.off() # Save it 
                
        } else {
                
                ##################################
                ######    NegaBin models   #######
                ##################################
          
                setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/binomial/",deparse(substitute(x)),"/"))
                calibrate.data=sample(1:nrow(dfplot2[!is.na(dfplot2$sp.mort.bin),]), Yportion*nrow(dfplot2[!is.na(dfplot2$sp.mort.bin),])) # 66% of our data 
                df1=dfplot2[!is.na(dfplot2$sp.mort.bin),][calibrate.data,]
                df2=dfplot2[!is.na(dfplot2$sp.mort.bin),][-calibrate.data,] 
                
                M1 = update(x, data=df1)
                M2 = predict (M1, newdata=df2,re.form=NA,type="response")
                
                ### Cor Test ###
                df2$names <- rownames(df2)
                M2 <- as.data.frame(M2)
                M2$names <-rownames(M2)
                df2 <- semi_join(df2,M2,by="names")
                D <- cor.test(df2$sp.mort.bin, M2$V1, method = ("pearson")) # To record somewhere
                capture.output(print(D), file=paste0(deparse(substitute(x)),"_CrossValid.Table.txt")) # Output as a latex wrapped in a txt file
                
                if (AllInOne == T){jpeg(file=paste0(deparse(substitute(x)),"Diagnostic.jpeg"),width=710);par(mfrow = c(2,2))
                } else {jpeg(file=paste0(deparse(substitute(x)),"_Response_VS_predicted.jpeg"),width=710)}
                plot(df2$sp.mort.bin, M2$V1,xlab="Fitted values", ylab="Predicted", main="Predicted vs. fitted",
                     cex.main=1.5, cex.lab =1.3) # To save 
                if (AllInOne == F){dev.off()} # To save 
                
                ### Residuals ###
                ResM <- residuals(M1)
                FitM <- fitted(M1)
                PredM <- predict(M1)
                if (AllInOne == F) jpeg(file=paste0(deparse(substitute(x)),"_Residuals_VS_fitted.jpeg"),width=710) # Save it 
                par(mfrow=c(1,1))
                plot(ResM ~ FitM, xlab="Fitted values", ylab="Residuals", main="Residuals vs. fitted",
                     cex.main=1.5, cex.lab =1.3)
                abline(h=0)
                if (AllInOne == F){dev.off()} # Save it 
                
                # Save this plot 
                if (AllInOne == F){jpeg(file=paste0(deparse(substitute(x)),"_QQplot.jpeg"),width=710)} # Save it 
                qqnorm(ResM, cex.main=1.5, cex.lab =1.3)
                qqline(ResM)
                if (AllInOne == F){dev.off()} # Save it 
                
                # Save this one 
                if (AllInOne == F){jpeg(file=paste0(deparse(substitute(x)),"_BinnedPlot.jpeg"),width=710)} # Save it 
                binnedplot(PredM,ResM)
                dev.off()
        }
        
        # Plot residuals against all my explicative variable
        c1 <- rep(ResM,times=length(ExplVar))
        c2 <- gather(M1$data[,c(ExplVar)])
        test <- cbind(c1,c2[,c(2,1)])
        jpeg(file=paste0(deparse(substitute(x)),"Residuals_Scatterplots.jpeg"),width=710)
        par(mfrow=c(1,1))
        xyplot(test[,1]~test[,2] | test[,3],col=1,
               strip=function(bg="white", ...)
                       strip.default(bg="white", ...),
               scales=list(alternating=TRUE,
                           x=list(relation="free"),
                           y=list(relation="same")),
               xlab="Explanatory variables",
               ylab="Response variable",
               panel=function(x,y){
                       panel.grid(h=-1,v=2)
                       panel.points(x,y,col=1)
                       panel.loess(x,y,xol=1,lwd=2)})
        dev.off()
        
        jpeg(file=paste0(deparse(substitute(x)),"_VisualisationTEST.jpeg"),width=710) # Save it 
        par(mfrow=c(3,3))
        par(ask=F)
        plot(x)
        dev.off()
        #plot(x)
        #dev.print(file=paste0(deparse(substitute(x)),"_Visualisation2.jpeg"),device=jpeg,width=710) # Save it 
        #dev.off()
        #dev.print(file=paste0(deparse(substitute(x)),"_Visualisation.jpeg"),device=jpeg,width=710) # Save it 
        #par(ask=F)
        #par(mfrow = c(1,1))
        
        
}
