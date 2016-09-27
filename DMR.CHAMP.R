DMR.CHAMP<-function(){
  library("ChAMP")
  Rlt<-list()
  samplesheetFile=list.files(pattern="*csv")
  ProjectName=unlist(strsplit(samplesheetFile,"[.]"))[1]
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
  myLoad<-champ.load(directory = getwd(),filterDetP=F,filterBeads=F,QCimages = F,filterXY=F,filterSNPs=F,filterMultiHit=F)
  # champ.CNA(intensity = myLoad$intensity, pd = myLoad$pd, loadFile = FALSE, batchCorrect = TRUE,
  #           file = "intensity.txt", resultsDir = paste(getwd(), "resultsChamp", sep = "/"),
  #           sampleCNA=TRUE, plotSample=TRUE, filterXY = TRUE, groupFreqPlots=TRUE,freqThreshold=0.3,
  #           control=TRUE,controlGroup="Control")
  # turn off the plot if you run it on the server since of the problem of X11
  myNorm<-champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet, pd = myLoad$pd, mset = myLoad$mset,sampleSheet = samplesheet, resultsDir = paste(getwd(), "resultsChamp",sep = "/"), methValue = "B", fromIDAT = TRUE, norm = "BMIQ", fromFile = FALSE, betaFile,filter = FALSE, filterXY = FALSE, QCimages = F, plotBMIQ = F)
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
  Fileout1<-paste(ProjectName,"resultsChamp.tar.gz",sep=".")
  Fileout2<-paste(ProjectName,"result.txt",sep=".")
  
  ## Manhattan Plot
  ManhattanPlot(mylimma)
  
  ## PCA Analysis
  PCAPlot(t(myNorm$beta),myLoad$pd$Sample_Group,ProjectName,multifigure=T)
  
  ## Cluster Analysis
  pheno<-myLoad$pd[match(colnames(myNorm$beta),myLoad$pd$Sample_Name),]$PublishID
  ClusterAnalysisPlot(myNorm$beta,pheno,suffix=ProjectName)
  
  ## TSNE ANLAYIS
  phen<-myLoad$pd[match(colnames(myNorm$beta),myLoad$pd$Sample_Name),]$Sample_Group
  TSNEAnalysis(t(myNorm$beta),pheno,prefix="TwinSomie")
  
  ## Heatmap Plot to DMR regions
  phen<-myLoad$pd[match(colnames(myNorm$beta),myLoad$pd$Sample_Name),]$Sample_Group
  HeatMap(myNorm$beta,phen=phen,varselect=1000,plot="heatmap.pdf",cexRow = 0.01,cexCol = 1.2,Colv=T,Rowv=T)
  
  ## MDS analysis
  pheno<-myLoad$pd[match(colnames(myNorm$beta),myLoad$pd$Sample_Name),]$PublishID
  MDSPlot(t(myNorm$beta),phen=pheno,prefix="TwinSmoke")
  
  cmd=paste("tar czvf ",Fileout1," ./resultsChamp",sep="");  
  system(cmd)
  #  write.table(result2,file=Fileout2,col.names=NA,row.names=T,sep="\t",quote=F)
  Rlt$dmr=result2
  Rlt$myNorm=myNorm
  Rlt$mylimma=mylimma
  Rlt$myLoad=myLoad
  return(Rlt)
}
