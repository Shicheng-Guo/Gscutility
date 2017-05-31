########################################################################################
###   Title: Tissue specific hyper-methylation biomarkers
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   updata time: 11/9/2015
###   Strategy: 1) Plasma based Cancer biomaker require tissue specific hyper-methylation. 
###   Strategy: 2) Methylation Biomarker should be low and non-methylated in PBMC samples
########################################################################################

saminfo<-read.table("../PancancerMethSaminfo_March2016.txt",sep="\t")
colname<-read.table("header.txt",sep="\t",row.names=1,colClasses=c("character",rep("character",6440)),nrow=1)
data<-read.table("PancancerMethMatrixjb",sep="\t",colClasses=c("character",rep("numeric",6440)),nrow=2500,row.names=1)
group<-paste(saminfo[match(colname,saminfo[,1]),3],"_",saminfo[match(colname,saminfo[,1]),4],sep="")

#' adopted group specific index
#' matrix are require at least one non-NA values in each row
#' @param data.frame of a matrix with colnames as the variable of group name
#' @return a data frame including GSI, reference assigned and AMF for each reference
#' @export
gsi<-function (matrix){
  narow<-which(unlist(apply(matrix,1,function(x) all(is.na(x)))))
  if(length(narow)>1){
    matrix<-matrix[-narow,]
  }
  group = names(table(colnames(matrix)))
  index = colnames(matrix)
  GSI <- c()
  gmaxgroup <- c()
  refMean<-c()
  refMax<-c()
  for (i in 1:nrow(matrix)) {
    gsit <- 0    
    refmean<-tapply(as.numeric(matrix[i, ]), index, function(x) mean(x, na.rm = T))
    refmax<-refmean[which.max(refmean)]
    gmax <- names(refmax)
    for (j in 1:length(group)) {
      tmp <- (1 - 10^(mean(na.omit(as.numeric(matrix[i, which(index == group[j])])), na.rm = T))/10^(mean(na.omit(as.numeric(matrix[i,which(index == gmax)])))))/(length(group) - 1)
      gsit <- c(gsit, tmp)
    }
    if(sum(refmean>0.3,na.rm=T)>1){
      gmax<-gsub(" ","",paste(unique(c(gmax,names(refmean)[which(refmean>0.3)])),',',collapse =""))
    }
    gmaxgroup <- c(gmaxgroup, gmax)
    GSI <- c(GSI, sum(gsit, na.rm = T))
    refMean<-c(refMean,gsub(" ","",paste(round(refmean,5),',',collapse ="")))
    refMax<-c(refMax,refmax)    
  }  
  if(nrow(matrix)!=length(refmax)){
    print("input matrix can not be allowed to be all NA in some row! please check your matrix!")
  }
  rlt = data.frame(region = rownames(matrix), group = gmaxgroup, GSI = GSI, refMax=refMax,refMean=refMean)
  return(rlt)
}

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    data1<-data[-NaRAW,]
  }else{
    data1<-data;
  }
  data1
}   

library("impute")
data<-read.table(file="outputFile",sep="\t",head=T,colClasses=c("character",rep("numeric",6440)),nrows=5,row.names=1,comment.char ="",as.is=T,check.names=F)


colnames(data)
strsplit()
data<-matrix(rnorm(485579*6500,1,1),485579,6500)
library(ff)
ffdf1 <- read.table.ffdf(file ="outputFile",header = TRUE,sep = "\t",nrows=20)


