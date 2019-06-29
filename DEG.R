DGE<-function(rnaseqdata,pid="TCGA-BRCA"){
yphen<-subset(phen,(bin==1 | bin==11) & pid==pid)
input<-rnaseqdata[,match(yphen$file_name,colnames(rnaseqdata))]
dim(yphen)
dim(input)
colnames(input)<-yphen$phen4
Seq<-paste(yphen$pid,yphen$bin,sep="-")
y<-abs(as.numeric(as.factor(yphen$bin))-2)
rlt<-c()
for(i in 1:nrow(input)){
  fit<-summary(bayesglm(y~input[i,]))$coefficients[2,]
  rlt<-rbind(rlt,fit)
}
rownames(rlt)<-rownames(input)
write.table(rlt,file=paste("~/hpc/methylation/",pid,"-RNAseq-FPKM-UQ.ENSG.DEG.txt",sep=""),sep="\t",col.names = T,row.names=T,quote=F)
rownames(rlt)<-ENSG2Symbol(rownames(rlt))
write.table(rlt,file=paste("~/hpc/methylation/",pid,"-RNAseq-FPKM-UQ.Symbol.DEG.txt",sep=""),sep="\t",col.names = T,row.names=T,quote=F)
}
