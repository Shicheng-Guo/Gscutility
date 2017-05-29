#!/usr/bin/env Rscript
library("monod")
library("ggplot2")
library("reshape2")
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input MHL matrix", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

Zscore<-function(ccpr,npr){
  Zmax<-matrix(0,nrow=nrow(ccpr),ncol=ncol(ccpr))
  for(i in 1:ncol(ccpr)){
    for(k in 1:nrow(npr)){
      idx<-1:ncol(npr)
      zmp <- (ccpr[k,i] - mean(npr[k, idx]))/(sd(npr[k, idx])*sqrt((length(npr[k, idx])-1)/(length(npr[k, idx]))))
      Zmax[k,i]=zmp 
    }
  }
  rownames(Zmax)=rownames(ccpr)
  colnames(Zmax)=colnames(ccpr)
  Zmax
}

prediction<-function(data,bio,tt=0.6){
  input<-data[na.omit(match(rownames(bio),rownames(data))),]
  bio<-bio[na.omit(match(rownames(input),rownames(bio))),]
  testcounts<-apply(input,2,function(x) table(unlist(lapply(bio[match(rownames(input)[which(x>tt)],rownames(bio)),]$group,function(x) unlist(strsplit(x,","))))))
  backcounts<-table(unlist(lapply(bio[match(rownames(input),rownames(bio)),]$group,function(x) unlist(strsplit(x,",")))))
  prediction<-apply(apply(testcounts,2,function(x) x/backcounts),2,function(x) rownames(backcounts)[which.max(x)])
  return(prediction)
}

load("/media/Home_Raid1/shg047/NAS1/monod/hapinfo/gsirlt.RData")
gsirlt1<-subset(gsirlt,refMax>0.6 & GSI>0.5)
gsirlt1$group<-as.character(gsirlt1$group)
bio<-gsi2bio(gsirlt1)
data=read.table(opt$file,head=T,row.names=1,sep="\t")

rlt<-c()
for(i in seq(0.2,0.8,0.1)){
  rlt<-rbind(rlt,prediction(data,bio,tt=i))
}
rownames(rlt)<-paste("R",seq(0,0.8,0.1),sep="")
print(t(data.frame(rlt)))

