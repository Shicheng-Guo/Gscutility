#  this code is for differential methylation region detection

setwd("/home/sguo/methylation/")
library("stringr")

PairWilPValue<-function(data,x1,x2){
  output<-matrix(NA,dim(data)[1],5)   # set output matrix ()
  for(i in 1:dim(data)[1]){
    out<-data.frame()
    if(all(! any(all(is.na(data[i,x1])),all(is.na(data[i,x2]))),sum(is.na(data[i,]))<0.5*length(data[i,]))){ 
      tmp1<-try(wilcox.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]),paired=T, na.action=na.omit))
      output[i,1]<-tmp1$p.value
      output[i,2]<-mean(as.numeric(data[i,x1]))-mean(as.numeric(data[i,x2]))
      output[i,3]<-"wilcox"
      output[i,4]<-mean(as.numeric(data[i,x1]))
      output[i,5]<-mean(as.numeric(data[i,x2]))
#      print(i)
    }
  }
  out<-cbind(rownames(data),output)
  out
}

result<-c()
P<-c()
sta<-c()

cancertype<-"LUAD"

cancertype<-rev(c("KIRC","BRCA","THCA","HNSC","PRAD","LIHC","KIRP","LUSC","COAD","UCEC","LUAD"))


for (cancer in rev(c("KIRC","BRCA","THCA","HNSC","PRAD","LIHC","KIRP","LUSC","COAD","UCEC","LUAD"))){
  pattern=paste("jhu-usc.edu_",cancer,".*",sep="")
  print (pattern)
  file=list.files(pattern=pattern)
  idv<-unique(as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*")))
  pairidv<-c()
  for (i in 1:length(idv)){
    t1<-paste(idv[i],"-01",sep="") 
    t2<-paste(idv[i],"-11",sep="")
    if(all(any(grepl(t1,file)),any(grepl(t2,file)))){
      pairidv<-c(pairidv,t1,t2)
    }
  }
  
  pairfile<-file[sapply(pairidv,function(x){x<-grep(x,file);x[length(x)]})]
  data<-c()
  for(i in 1:length(pairfile)){
    tmp<-read.table(pairfile[i],head=T,sep="\t",as.is=F,skip=1)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
    data<-cbind(data,tmp[,2])
    print(c(i,cancer))
  }
  rownames(data)=tmp[,1]
  colnames(data)=pairidv
  RData<-paste(cancer,"meth450.pair.RData",sep="")
  save(data,file=RData)
  # feature selection and classification
  # differential methylation calculation
  data<-data+matrix(rnorm(length(data),0.0001,0.0001),dim(data)[1],dim(data)[2])   # row=gene, col=inv
  type<-substr(colnames(data),14,15)
  x1<-which(type==names(table(type))[1])   # type 1, cancer or sensitive
  x2<-which(type==names(table(type))[2])   # type 2, normal or resistant
  out<-PairWilPValue(data,x1,x2)
  colnames(out)<-paste(cancer,c("GeneName","Pvalue","Statistic","test","mean1","mean2"),sep=".")
  P<-cbind(P,out[,2])
  sta<-cbind(sta,out[,3])
  result<-cbind(result,out)
}

Pvalue=data.frame(tmp[,3],P)
colnames(Pvalue)=c("GeneSymbol",paste(cancertype,c("Pvalue"),sep="."))
result=data.frame(tmp[,3],result)
colnames(result)=c("GeneSymbol",paste(cancertype,c("GeneName","Pvalue","Statistic","test","mean1","mean2"),sep="."))

write.table(Pvalue, file="PanCancerPvalue.txt",sep="\t",col.names=NA,row.names=T)
write.table(result, file="PanCancerFullResult.txt",sep="\t",col.names=NA,row.names=T)


# The pattern of P


