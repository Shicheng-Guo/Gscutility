library("GEOquery")
# GSE52826 <- getGEO("GSE52826",destdir="/home/sguo/monod/data/geo")
# save(GSE52826, file="GSE53045_matrix.Rdata")
load("GSE53045_matrix.Rdata")
data <- as.data.frame(exprs(GSE52826[[1]]))
phen <- pData(phenoData(GSE52826[[1]]))
# type1<-phen[match(colnames(data),rownames(phen)),]$characteristics_ch1
# PCAPlot(t(na.omit(data)),type1,output="esca.pca.GSE52826.pdf",legend.cex=0.5,multifigure=T)
# type2<-phen[match(colnames(data),rownames(phen)),]$characteristics_ch1.1
library("ChAMP")
newdata<-na.omit(data)
dim(newdata) # 479037/485577=98.6%
pd<-c()
pd$Sample_Name=as.character(rownames(phen))
pd$Sentrix_ID=1:nrow(phen)
pd$Sentrix_Position=1:nrow(phen)
pd$Sample_Group=as.character(unlist(lapply(phen$characteristics_ch1.1,function(x) unlist(strsplit(as.character(x),": "))[2])))
pd<-data.frame(pd)
myNorm=na.omit(data)
colnames(myNorm)=pd$Sample_Name
beta.norm=myNorm
compare.group = c("tumor","normal")
mylimma=champ.MVP(beta.norm =myNorm, pd =pd, adjPVal = 0.05, adjust.method = "BH",compare.group = c("tumor","normal"), resultsDir = paste(getwd(), "resultsChamp", sep = "/"),bedFile = TRUE)
library("doParallel")
champ.DMR(betaNorm=data.matrix(myNorm),design=pd$Sample_Group,maxGap=300,
          cutoff=0.2,minProbes=3,smooth=TRUE,smoothFunction=loessByCluster,
          useWeights=FALSE,permutations=NULL,B=100,pickCutoff=FALSE,
          pickCutoffQ=0.99,nullMethod="bootstrap",verbose=TRUE,cores=5,arraytype="450K",
          method = "Bumphunter",resultsFile=mylimma,meanLassoRadius=375,minSigProbesLasso=3,minDmrSep=1000,
          minDmrSize=50,adjPvalProbe=0.05,adjPvalDmr=0.05,pData=pd)
