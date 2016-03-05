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



gsi<-function(data){
  group=names(table(colnames(data)))
  index=colnames(data)
  gsi<-c()
  gmaxgroup<-c()
  for(i in 1:nrow(data)){
    gsit<-0
    gmax<-names(which.max(tapply(as.numeric(data[i,]),index,mean)))
    for(j in 1:length(group)){
      tmp<-(1-10^(mean(data[i,][which(index==group[j])]))/10^(mean(data[i,][which(index==gmax)])))/(length(group)-1)
      gsit<-gsit+tmp
    }
    gmaxgroup<-c(gmaxgroup,gmax)
    gsi<-c(gsi,gsit)
    print(c(gmax,gsit))
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi)
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


