# Alex on the 23/05/2018
# Script to save the output of my models with package spaMM
# It also requires the function Rsquared

SavingInfo = "Saving \n A R function to save your model and various output  \n Require the knitr & the spaMM packages. RandomFields might be needed if there is spatial autocorrelation \n \n Input data: Your model (HLfit) \n 
create and save four elements: \n (1) A directory named after your model \n (2) A csv file with the output results of the Rsquared.AC function \n (3) A csv file with the estimated coefs, standards deviations and significance levels \n (4) A txt file with these information as LaTex format \n (5) The pattern of the spatial autocorrelation is also plotted "
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

library(knitr)
library(spaMM)
library(fields)
library(piecewiseSEM)

Saving <- function(x){
  dir.create(paste0(Dir,"/",deparse(substitute(x)),"/"))
  save(x, file = paste0(Dir,"/",deparse(substitute(x)),"/",deparse(substitute(x)),".rda")) # save the model as an RDA file 
  if (x$family$family=="binomial"){Namecat <- "Binomial"} else Namecat <- "Negbin" # Identidy the family
  y <- data.frame(matrix(unlist(rsquared.AC(x)), nrow=1,byrow=T,dimnames = list(c(NULL),c("family","link","method","Marginal","Conditional","Lik","AICm","Craw1","ICCadj1", "PCVran", "PCVObs","ICCraw2","ICCadj2","PCV1","PCV2"))))
  if (length(list.files(path=paste0(Dir,"/"),pattern=paste0("Models_",Namecat,CODE,"_",seuil,".csv")))==0){
    write.table(y,paste0(Dir,"/Models_",Namecat,CODE,"_",seuil,".csv"),append=T, quote = FALSE, sep=";",row.names = deparse(substitute(x)),col.names=NA) # If file does not exist
  }else write.table(y, paste0(Dir,"/Models_",Namecat,CODE,"_",seuil,".csv"),append=T, quote = FALSE, sep=";",row.names = deparse(substitute(x)),col.names=F)  # If it does exist                             
  # The rsquared output are stored at the end of each other in a huge Csv file
  
  A <- as.data.frame(summary(x)[2]) # All intercepts fixed
  
  if (length(grep("Matern",x[["predictor"]],fixed=T))>=1)    # If there is a spatial effect
  {B <- cbind("",t(t(as.numeric(unclass(x[["lambda"]])))),"") # Extract the lambda values in a table 
  colnames(B) <- colnames(A)
  rownames(B) <- sub(" .","", names(x[["lambda"]]), ignore.case = FALSE,fixed = T) # Name them correctly
  B <- rbind(B,B) # Duplicate this table twice 
  B[1:(nrow(B)/2),1] <- "Log Lambda" # First half of the first column are loglambda
  B[(1+(nrow(B)/2)):nrow(B),1] <- "Lambda" # The second half of the first column is lambda
  B[1:(nrow(B)/2),2] <- log(as.numeric(B[(1+(nrow(B)/2)):nrow(B),2])) # The second half of the second column is then equal to the exponential of the first half
  D <- cbind("Estimates",t(t(unlist(x$ranFix[[1]]))),"") # Extract the estimated parameters of the spatial effect
  rownames(D) <- sub("2.","",rownames(D), ignore.case = FALSE,fixed = T)
  B <- rbind(B,D) # Combine all together
  par(mfrow=c(1,1))
  d<- seq(0,40,length.out=80)
  y<- Matern(d, range=as.numeric(D[2,2]), smoothness=as.numeric(D[1,2])) # Plot my spatial effects
  matplot( d, y, type="l", lty=1, lwd=2,main=c(bquote(paste("Estimated parameters : ",rho," = ",.(round(as.numeric(D[2,2]),3))," and ",nu," = ",.(round(as.numeric(D[1,2]),3))))),ylab="Response covariance",xlab="Distance between plots")
  dev.print(file=paste0(Dir,"/",deparse(substitute(x)),"/",deparse(substitute(x)),"Autocor.jpeg"),device=jpeg,width=710) # Save it 
  
  } else { # If there is no spatial effect 
    
    B <- as.data.frame(summary(x)[3])[,2:4] # Extract the random effects 
    colnames(B) <- colnames(A)
    B <- rbind(B,B) # Duplicate this table twice
    B[,1] <- as.character(B[,1]) # Make the first column as a character
    B[1:(nrow(B)/2),1] <- "Log Lambda" # First half of the first column are loglambda
    B[(1+(nrow(B)/2)):nrow(B),1] <- "Lambda" # The second half of the first column is lambda
    B[(1+(nrow(B)/2)):nrow(B),2] <- exp(B[1:(nrow(B)/2),2]) # The first half of the second column is then equal to the log of the second half
    B[(1+(nrow(B)/2)):nrow(B),3] <- "" # The second half of the third column is empty
  }
  
  C <- rbind(A,B) # All together (random + fixed effects)
  colnames(C) <- c("Estimate","Cond. SE","t-value")
  C[,"Signif"] <- ""
  C[C$`t-value`!="" & (as.numeric(C$`t-value`)>=2 | as.numeric(C$`t-value`)<=-2),"Signif"] <- "*" # Add significiance 
  write.table(C,paste0(Dir,"/",deparse(substitute(x)),"/",deparse(substitute(x)),"_Table.csv"),sep=";",col.names = NA, row.names = T) # Variable outputs in a table.csv
  D <- kable(C,digits = 3,"latex",align=c("rrrc"))
  capture.output(print(D), file=paste0(Dir,"/",deparse(substitute(x)),"/",deparse(substitute(x)),"_Tex.Table.txt")) # Output as a latex wrapped in a txt file
}

