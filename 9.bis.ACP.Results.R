Allcode <- c("QUEILE","PINPINA","PINSYL","PICABI","FAGSYL","PINHAL","QUEROB","PINNIG","QUEPET","CASSAT","ABIALB","QUEPYR","FRAEXC","PINPIN",
        "QUESUB","BETPEN","POPTRE","POPNIG")
ACP.all <- data.frame(matrix(nrow=21,ncol=3*18,NA))
names1 <- c("poids.tot","poids","poids.rank")
allnames <- expand.grid(names1,Allcode)
Tablenames <- paste0(allnames[,2],allnames[,1])
colnames(ACP.all) <- Tablenames

ACP.sum <- data.frame(matrix(nrow=18,ncol=6,NA))
colnames(ACP.sum) <- c("Species","Min Variable","Min %","Max Variable", "Max %","Mean %")

i = 1
for (code in Allcode){
acp10000 <- readRDS(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",code,"/PCA/acp10000.rds"))
ACP.all[,paste0(code,"poids.tot")] <- apply(acp10000$co,1,function(x) sum(abs(x)*c(acp10000$eig[1:2]/sum(acp10000$eig)*100)))
#poids2 <- poids[order(poids,decreasing=T)]
ACP.all[,paste0(code,"poids")] <- ACP.all[,paste0(code,"poids.tot")]/sum(ACP.all[,paste0(code,"poids.tot")])*100
ACP.all[,paste0(code,"poids.rank")] <- rank(-ACP.all[,paste0(code,"poids.tot")])

ACP.sum[i,1] <- code
ACP.sum[i,2] <- rownames(acp10000$co)[which(ACP.all[,paste0(code,"poids")]==min(ACP.all[,paste0(code,"poids")]))]
ACP.sum[i,3] <- min(ACP.all[,paste0(code,"poids")])
ACP.sum[i,4] <- rownames(acp10000$co)[which(ACP.all[,paste0(code,"poids")]==max(ACP.all[,paste0(code,"poids")]))]
ACP.sum[i,5] <- max(ACP.all[,paste0(code,"poids")])
ACP.sum[i,6] <- mean(ACP.all[,paste0(code,"poids")])
i = i+1
}

ACP.all[,"Variable"] <- rownames(acp10000$co)