# this code is for differential methylation sites based on Pan-cancer methylation data
# input data is PanmethylationData.Rdata (PancancerMethMatrix_March2016.RData)
# tradtional t-test was applied, it should be replaced by other in the comming days
# DMS in multiple cancer is quite interesting in current stage

setwd("/media/Home_Raid1/shg047/NAS3/HM450/TCGA/Data/split")
library("stringr")
javascript:void(0);
load("PancancerMethMatrix_March2016.RData")

PairTtestPValue<-function(data,x1,x2){
  output<-matrix(NA,dim(data)[1],7)   # set output matrix ()
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]),paired=T, na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output[i,2]<-mean(as.numeric(data[i,x1]))-mean(as.numeric(data[i,x2]))
      output[i,3]<-"ttest"
      output[i,4]<-mean(as.numeric(data[i,x1]))
      output[i,5]<-mean(as.numeric(data[i,x2]))
      output[i,6]<-sd(data[i,x1])
      output[i,7]<-sd(data[i,x2])
      print(i)
    }
  }
  out<-cbind(rownames(data),output)
  out
}

colnames(data)<-unlist(lapply(colnames(data),function(x) gsub("[.]","-",x)))
cancertype<-unique(unlist(lapply(colnames(data),function(x) unlist(strsplit(x,"_|-Human"))[2])))

result<-c()
P<-c()
sta<-c()
CANCER<-c()
for (cancer in cancertype){
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
  if(length(pairidv)>3){
  pairfile<-file[sapply(pairidv,function(x){x<-grep(x,file);x[length(x)]})]
  newdata<-data[,match(pairfile,colnames(data))]
  colnames(newdata)=pairidv
  # feature selection and classification
  # differential methylation calculation
  newdata<-newdata+matrix(rnorm(length(newdata),0.0001,0.0001),dim(newdata)[1],dim(newdata)[2])   # row=gene, col=inv
  type<-substr(colnames(newdata),14,15)
  x1<-which(type==names(table(type))[1])   # type 1, cancer or sensitive
  x2<-which(type==names(table(type))[2])   # type 2, normal or resistant
  out<-PairTtestPValue(data,x1,x2)
  colnames(out)<-paste(cancer,c("GeneName","Pvalue","Statistic","test","mean1","mean2","SD1","SD2"),sep=".")
  CANCER<-c(CANCER,cancer)
  P<-cbind(P,out[,2:3])
  sta<-cbind(sta,out[,3])
  result<-cbind(result,out)
  colnames(P)<-paste(rep(CANCER,each=2),c("P","Delta"),sep="-") 
  }
  outputfile=paste("TCGA-",cancer,".DMS.txt",sep="")
  write.table(out,file=outputfile,col.names=T,row.names=F,sep="\t",quote=F)
}

P<-data.frame(out[,1],P)
write.table(P, file="PanCancerPvalue.txt",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(result, file="PanCancerFullResult.txt",sep="\t",col.names=NA,row.names=T,quote=F)




