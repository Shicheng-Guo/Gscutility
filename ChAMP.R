# run the script in the fold of idat of HM450
# update in April 2016

library("ChAMP")
library("outliers")
setwd("/home/sguo/dyh/idat")
myLoad<-champ.load(directory = getwd(),filterBeads=TRUE,QCimages = F)
save(myLoad,file="myLoad.RData")
# turn off the plot if you run it on the server since of the problem of X11
myNorm<-champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet, pd = myLoad$pd, mset = myLoad$mset,sampleSheet = "sampleSheet.txt", resultsDir = paste(getwd(), "resultsChamp",sep = "/"), methValue = "B", fromIDAT = TRUE, norm = "BMIQ", fromFile = FALSE, betaFile,filter = TRUE, filterXY = TRUE, QCimages = F, plotBMIQ = F)
save(myNorm,file="myNorm.RData")
TMR<-champ.TrueMethyl(beta.norm = myNorm$beta, pd = myLoad$pd, adjPVal = 0.5, adjust.method = "BH",compare.group = c("C", "T"), resultsDir = paste(getwd(), "resultsChamp", sep = "/"),bedFile = TRUE)



# 2016-May-update
###################################################################################################################
############################ Illumina Methylation Microarray Analysis Pipeline ####################################
###################################################################################################################

#================================================================================================================
#=============================================  Input: Beta-Matrix  =============================================
#================================================================================================================

library("ChAMP")
data<-list()
beta=data.frame(abs(matrix(rnorm(10*6,0.1,0.2),10,5)),abs(matrix(rnorm(10*6,0.8,0.2),10,5)))
colnames(beta)=c(paste("T",1:5,sep=""),paste("C",1:5,sep=""))
rownames(beta)=c("cg00114008","cg00003287","cg00026222","cg00065215","cg00038675","cg00083652","cg00036263","cg00046018","cg00052626","cg00046625")
pd<-data.frame(matrix(NA,ncol(beta),10))
colnames(pd)=c("Sample_Name","Sample_Plate","Sample_Group","Pool_ID","Project","Sample_Well","Array","Slide","Basename","filenames")

head(beta)
pd$Sample_Name=colnames(beta)
pd$Sample_Group=c(rep("T",5),rep("C",5))
data$beta=beta
data$pd=pd
pd

# DMR analysis
result1=champ.MVP(beta.norm = data$beta, pd = data.frame(data$pd), adjPVal = 0.05, adjust.method = "BH",
          compare.group = c("C", "T"), resultsDir =getwd(),bedFile = TRUE)
write.table(result1,file="",sep="\t",col.names=NA,row.names=T,quote=F)

Result2<-champ.lasso(beta.norm = data$beta, pd = data.frame(data$pd), filterXY = TRUE, image = TRUE, mafPol.lower = 0,
            mafPol.upper = 0.05, popPol = "eur", lassoStyle = "max", lassoRadius = 200,
            minSigProbesLasso = 3, minDmrSep = 1000, minDmrSize = 0, adjPVal = 0.05,
            adjust.method = "BH", resultsDir = getwd(),
            bedFile = TRUE, DMRpval = 0.05, batchDone = FALSE)

#================================================================================================================
#=============================================  Input: idata ====================================================
#================================================================================================================

data(testDataSet)
myLoad=testDataSet
dim(myLoad$beta)
myNorm=champ.norm(norm="NONE")
rm(list=ls())





