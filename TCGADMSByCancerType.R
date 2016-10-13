##########################################################################################################################
###################################### Introduction and Usage  #########################################################
##########################################################################################################################

#  this code is for differential methylation region detection
#  Usage: Rscript PanMethBiomark.R "BLCA"
#  Date: 10/12/2016
#########################################################################################################################
###################################### Function and Package Load #########################################################
##########################################################################################################################

PairTtestPValue<-function(data,x1,x2){
  data<-data.matrix(data)
  output<-matrix(NA,dim(data)[1],6)   # set output matrix ()
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(t.test((data[i,x1]),(data[i,x2]),paired=T, na.action=na.omit))
      output[i,1]<-format(tmp1$p.value, scientific=TRUE)
      output[i,2]<-round(mean((data[i,x1]))-mean((data[i,x2])),3)
      output[i,3]<-round(mean((data[i,x1])),3)
      output[i,4]<-round(mean((data[i,x2])),3)
      output[i,5]<-round(sd(data[i,x1]),3)
      output[i,6]<-round(sd(data[i,x2]),3)
      print(i)
    }
  }
  out<-cbind(rownames(data),output)
  out
}


##########################################################################################################################
###################################### Working Pipeline #########################################################
##########################################################################################################################
args <- commandArgs(trailingOnly = TRUE)
print(args)

setwd("/media/Home_Raid1/shg047/NAS3/HM450/TCGA/Data/split")
library("stringr")
load("PancancerMethMatrix_March2016.RData")
#load("PancancerMethMatrix_March2016.Test.RData")
colnames(data)<-unlist(lapply(colnames(data),function(x) gsub("[.]","-",x)))
# cancertype<-unique(unlist(lapply(colnames(data),function(x) unlist(strsplit(x,"_|-Human"))[2])))

cancer<-args[1]
CANCER<-c()
sta<-c()
result<-c()
P<-c()

# Identify Paired Tumor-Adjacent Samples(01/11) 
file=colnames(data)[grep(cancer,colnames(data))]
idv<-unique(as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*")))
pairidv<-c()
for (i in 1:length(idv)){
  t1<-paste(idv[i],"-01",sep="") 
  t2<-paste(idv[i],"-11",sep="")
  if(all(any(grepl(t1,file)),any(grepl(t2,file)))){
    pairidv<-c(pairidv,t1,t2)
  }
}

# Pair-wise DMS test
if(length(pairidv)>3){
  pairfile<-file[sapply(pairidv,function(x){x<-grep(x,file);x[length(x)]})]
  newdata<-data[,match(pairfile,colnames(data))]
  colnames(newdata)=pairidv
  newdata<-newdata+matrix(rnorm(length(newdata),0.0001,0.0001),dim(newdata)[1],dim(newdata)[2])   # row=gene, col=inv
  type<-substr(colnames(newdata),14,15)
  x1<-which(type==names(table(type))[1])   # type 1, cancer or sensitive
  x2<-which(type==names(table(type))[2])   # type 2, normal or resistant
  Rlt<-PairTtestPValue(newdata,x1,x2)
}

Rlt
FDR<-p.adjust(Rlt[,2],method="fdr")
Rlt<-data.frame(Rlt,format(FDR, scientific=TRUE))
colnames(Rlt)<-paste(cancer,c("GeneName","Pvalue","Statistic","mean1","mean2","SD1","SD2","FDR"),sep=".")
outputfile=paste("TCGA-",cancer,".DMS.txt",sep="")
write.table(Rlt,file=outputfile,col.names=T,row.names=F,sep="\t",quote=F)







