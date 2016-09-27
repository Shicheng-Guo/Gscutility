# use latest champ in lastest R version. 
# Sep 26th 2016
# We obtain 12 pairs twins mh450K array (smoking vs non-smoking)

DMR.CHAMP<-function(){
  library("ChAMP")
  Rlt<-list()
  samplesheetFile=list.files(pattern="*csv")
  sampleID=unlist(strsplit(samplesheetFile,"[.]"))[1]
  # be sure to remove header of the samplesheet
  sampleinfo<-read.table(samplesheetFile,head=T,sep="\t",as.is=T)
  sampleinfo
  system("mkdir resultsChamp")
  # champ.process(fromIDAT = TRUE, fromFile = FALSE, directory = getwd(), resultsDir =paste(getwd(), "resultsChamp", sep = "/"), methValue = "B", filterDetP = TRUE,
  #              detPcut = 0.01, filterXY = TRUE, removeDetP = 0, filterBeads = TRUE, beadCutoff =
  #                0.05, filterNoCG = FALSE, QCimages = TRUE, batchCorrect = FALSE, runSVD =
  #                TRUE, studyInfo = FALSE, infoFactor = c(), norm = "BMIQ", adjust.method = "BH",
  #              adjPVal = 0.05, runDMR = TRUE, runCNA = FALSE, plotBMIQ = FALSE, DMRpval = 0.05,
  #              sampleCNA=FALSE,plotSample = TRUE,groupFreqPlots=TRUE,freqThreshold=0.3, bedFile
  #              = TRUE, methProfile = TRUE, controlProfile = FALSE)
  myLoad<-champ.load(directory = getwd(),filterBeads=TRUE,QCimages = F)
  # champ.CNA(intensity = myLoad$intensity, pd = myLoad$pd, loadFile = FALSE, batchCorrect = TRUE,
  #           file = "intensity.txt", resultsDir = paste(getwd(), "resultsChamp", sep = "/"),
  #           sampleCNA=TRUE, plotSample=TRUE, filterXY = TRUE, groupFreqPlots=TRUE,freqThreshold=0.3,
  #           control=TRUE,controlGroup="Control")
  # turn off the plot if you run it on the server since of the problem of X11
  myNorm<-champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet, pd = myLoad$pd, mset = myLoad$mset,sampleSheet = samplesheet, resultsDir = paste(getwd(), "resultsChamp",sep = "/"), methValue = "B", fromIDAT = TRUE, norm = "BMIQ", fromFile = FALSE, betaFile,filter = TRUE, filterXY = TRUE, QCimages = F, plotBMIQ = F)
  # identify MVP and it would be used in DMR calling
  mylimma=champ.MVP(beta.norm = myNorm$beta, pd = data.frame(myLoad$pd), adjPVal = 0.5, adjust.method = "BH",compare.group = c("R","S"), resultsDir = paste(getwd(), "resultsChamp", sep = "/"),bedFile = T)
  # Do DMR calling and save it to result2
  library("doParallel")
  result2=champ.DMR(betaNorm=myNorm$beta,design=myLoad$pd$Sample_Group,maxGap=300,
                    cutoff=0.1,minProbes=3,smooth=TRUE,smoothFunction=loessByCluster,
                    useWeights=FALSE,permutations=NULL,B=10,pickCutoff=FALSE,
                    pickCutoffQ=0.99,nullMethod="bootstrap",verbose=TRUE,cores=3,arraytype="450K",
                    method = "Bumphunter",resultsFile=mylimma,meanLassoRadius=375,minSigProbesLasso=3,minDmrSep=1000,
                    minDmrSize=50,adjPvalProbe=0.5,adjPvalDmr=0.5,pData=myLoad$pd)
  Fileout1<-paste(sampleID,"resultsChamp.tar.gz",sep=".")
  Fileout2<-paste(sampleID,"result.txt",sep=".")
  cmd=paste("tar czvf ",Fileout1," ./resultsChamp",sep="");  
  system(cmd)
  #  write.table(result2,file=Fileout2,col.names=NA,row.names=T,sep="\t",quote=F)
  Rlt$dmr=result2
  Rlt$myNorm=myNorm
  Rlt$mylimma=mylimma
  Rlt$myLoad=myLoad
  return(Rlt)
}

# prepare csv samplesheet and each samplesheet will implement the prevous DMR scrpt
# if you creat multiple csv to conduct multiple different DMR analysis but be sure just keep one samplesheet within idat folder. 
# if you have multiple csv samplesheet within idat or sub-idat folder, there will be error for the script. 
# For twins smoking dataset, we only have one comparsion way, therefore, only one csv samplesheet. 
# copy all the idat and corresponding cvs samplesheet into a new folder (here, idat)

setwd("/media/NAS3_volume2/shg047/methage/twin/idat")
Rlt<-DMR.CHAMP()
image=paste(colnames(R)[i],"image.RData",sep=".")
save.image(file=image)
