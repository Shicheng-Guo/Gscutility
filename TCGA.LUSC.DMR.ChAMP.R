#######################################################################################################################
###   Title : DMR and DVR analysis to lusc dataset in TCGA Project with CHAMP packages
###   Author: Shicheng Guo, Ph.D. Email: Shicheng.Guo@hotmail.com 
###   Time :  Sep/11/2015 
###   Section 1.  Download DNA methylation and RNA-seqV2 dataset, clinical information of lusc and lusc
###   Section 2.  Annotation and databased pre-loading 
###   Section 3.  CHAMP Pipeline Analysis (MVP and DMR)
###   Section 4.  Summary
#######################################################################################################################

#######################################################################################################################
#######################################  library loading ##############################################
#######################################################################################################################

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


DGEtest<-function(data,phen,compare.group=c("01","11")){
  output<-matrix(NA,dim(data)[1],5)   # set output matrix ()
  x1<-which(phen==compare.group[1])   # type 1, cancer or sensitive
  x2<-which(phen==compare.group[2])   # type 2, normal or resistant
  a.shapiro.pvalue<-c()
  b.shapiro.pvalue<-c()
  delta<-mean.A<-mean.B<-tmp.ttest<-tmp.wilcox<-c()
  Pvalue.ttest<-Pvalue.wilcox<-c()
  for(i in 1:nrow(data)){
    tmp.a.shapiro.pvalue<-try(shapiro.test(as.numeric(data[i,x1]))$p.value)
    tmp.b.shapiro.pvalue<-try(shapiro.test(as.numeric(data[i,x2]))$p.value)
    a.shapiro.pvalue<-c(a.shapiro.pvalue,tmp.a.shapiro.pvalue)
    b.shapiro.pvalue<-c(b.shapiro.pvalue,tmp.b.shapiro.pvalue)
    delta<-c(delta,mean(as.numeric(data[i,x1]))-mean(as.numeric(data[i,x2])))
    mean.A<-c(mean.A,mean(as.numeric(data[i,x1])))
    mean.B<-c(mean.B,mean(as.numeric(data[i,x2])))
    tmp.ttest<-try(t.test(as.numeric(data[i,x1]),as.numeric(data[i,x2],paired=F)))
    Pvalue.ttest<-c(Pvalue.ttest,tmp.ttest$p.value)
    tmp.wilcox<-try(wilcox.test(as.numeric(data[i,x1]),as.numeric(data[i,x2]),paired=F))
    Pvalue.wilcox<-c(Pvalue.wilcox,tmp.wilcox$p.value)
  }
  p.adjust.BH<-p.adjust(Pvalue.ttest,"BH")
  output<-data.frame(mean.A,mean.B,delta,Pvalue.ttest,p.adjust.BH,Pvalue.wilcox,a.shapiro.pvalue,b.shapiro.pvalue)
  rownames(output)<-rownames(data)
  return(output)
}

library("impute")
library("stringr")
library("ChAMP")

#######################################################################################################################
#######################################  lusc DNA methylation Collection ##############################################
#######################################################################################################################
# primary directory
priDir="/home/sguo/tcga/lusc/"                                                        # main directory
methDir=paste(priDir,"DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/",sep="")  # Methylation Dataset
methDataSave1="TCGA.Methy450.lusc.RData"                                              # save methylation RData
methDataSave2="TCGA.Methy450.lusc.RM.Impute.RData"                                    # save imputate methylation RData
setwd(methDir)

file=list.files(pattern="jhu-usc.edu_*")
idv<-unique(as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*")))
idv1<-as.array(str_extract(file,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*"))
pairidv<-c()
for (i in 1:length(idv)){
  t1<-paste(idv[i],"-01",sep="") 
  t2<-paste(idv[i],"-11",sep="")
  if(all(any(grepl(t1,file)),any(grepl(t2,file)))){
    pairidv<-c(pairidv,t1,t2)
  }
}
allfile<-file[sapply(idv1,function(x){x<-grep(x,file);x[length(x)]})]
data1<-c()
for(i in 1:length(allfile)){
  tmp<-read.table(allfile[i],head=T,sep="\t",as.is=F,skip=1)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  data1<-cbind(data1,tmp[,2])
}
map<-tmp[,c(1,3,4,5)]
rownames(data1)=tmp[,1]
colnames(data1)=idv1
methdata<-list()
methdata$matrix<-data1
methdata$map<-map
save(methdata,file=paste(priDir,methDataSave1,sep=""))

methdata=RawNARemove(data1)
newmethdata<-impute.knn(methdata)$data
rownames(newmethdata)=rownames(methdata)
colnames(newmethdata)=colnames(methdata)
newmethdata<-newmethdata+matrix(rnorm(length(newmethdata),0.0001,0.0001),dim(newmethdata)[1],dim(newmethdata)[2])   # row=gene, col=inv
methdata<-list()
methdata$matrix<-newmethdata
methdata$map<-map[match(rownames(newmethdata),map[,1]),]
save(methdata,file=paste(priDir,methDataSave2,sep=""))

# summary
print("Methylation File reading and combinding completed!")
print("Methylation File reading and combinding is time-consuming, therefore, the RData has been stored when the data were read in the first time!")
print(paste("Totally,",length(idv1),"methylation files were enrolled from TCGA project. and the sample size is as the following:",sep=" "))
print(table(substr(idv1,14,15)))
print(paste("Among these samples,",length(pairidv),"samples are paired samples",sep=" "))

#######################################################################################################################
#######################################  lusc DNA methylation Analysis ##############################################
#######################################################################################################################
priDir="/home/sguo/tcga/lusc/"
methDataSave2="TCGA.Methy450.lusc.RM.Impute.RData"
setwd(priDir)
load(methDataSave2)
library("ChAMP")

beta=methdata$matrix
beta=methdata$matrix[,c(grep("01",substr(colnames(beta),14,15)),grep("11",substr(colnames(beta),14,15)))]
pd<-data.frame(matrix(NA,ncol(beta),ncol(beta)))
colnames(pd)=c("Sample_Name","Sample_Plate","Sample_Group","Pool_ID","Project","Sample_Well","Array","Slide","Basename","filenames")
pd$Sample_Name=colnames(beta)
pd$Sample_Group=substr(colnames(beta),14,15)
compare.group = c("01", "11")

data=list()
data$beta=beta
data$pd=pd
system("mkdir resultsChamp")
rlt1<-champ.MVP(beta.norm = data$beta, pd = data.frame(data$pd), adjPVal = 0.05, adjust.method = "BH",
                compare.group = compare.group, resultsDir = paste(getwd(), "resultsChamp", sep = "/"),
                bedFile = TRUE)
rlt2<-champ.lasso(beta.norm = data$beta, pd = data.frame(data$pd),filterXY = TRUE, image = TRUE, mafPol.lower = 0,
                  mafPol.upper = 0.05, popPol = "eur", lassoStyle = "max", lassoRadius = 2000,
                  minSigProbesLasso = 3, minDmrSep = 1000, minDmrSize = 0, adjPVal = 0.05,
                  adjust.method = "BH", resultsDir = paste(getwd(), "resultsChamp", sep = "/"),
                  bedFile = TRUE, DMRpval = 0.05, batchDone = FALSE)

write.table(rlt1,file="TCGA.Meth450.lusc.ChAMP.DMS.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(rlt2$dmrs,file="TCGA.Meth450.lusc.ChAMP.DMR.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(rlt2$dmr.probes,file="TCGA.Meth450.lusc.ChAMP.DMR.Probe.txt",col.names=T,row.names=F,quote=F,sep="\t")

#######################################################################################################################
#######################################  lusc Gene expression analysis ##############################################
#######################################################################################################################

setwd("/home/sguo/tcga/lusc/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3")
file=list.files(pattern="*rsem.genes.normalized_results")
sampleinfo=read.table("/home/sguo/tcga/lusc/file_manifest.txt",head=T,sep="\t")
expDataSave1="/home/sguo/tcga/lusc/TCGA.RNAseqV2.lusc.RData"    # save RNAseqV2 matrix
exptxtSave1="/home/sguo/tcga/lusc/TCGA.RNAseqV2.lusc.Diff.txt"  # save signficant differential gene ouput

sampleid=sampleinfo[match(file,sampleinfo[,ncol(sampleinfo)]),5]
expdata<-c()
for(i in 1:length(file)){
  tmp<-read.table(file[i],head=T,sep="\t",as.is=F)  # tmp<-read.table(file[i],head=T,sep="\t",as.is=F)
  expdata<-cbind(expdata,tmp[,2])
  print(i)
}
rownames(expdata)=tmp[,1]
colnames(expdata)=sampleid
save(expdata,file=expDataSave1)
type<-substr(sampleid,14,15)
table(type)
expdata<-expdata+matrix(rnorm(length(expdata),0.0001,0.0001),dim(expdata)[1],dim(expdata)[2])   # row=gene, col=inv
rlt1<-DGEtest(expdata,phen=type,compare.group=c("01","11"))
write.table(rlt1,file=exptxtSave1,col.names=T,row.names=T,quote=F,sep="\t")

#######################################################################################################################
################################# DNA methylation and Gene expression integrate ##############################
######################################################################################################################
setwd(priDir)

methy.data<-read.table("TCGA.Meth450.lusc.ChAMP.DMR.txt",head=T,sep="\t",as.is=T)
exp.data<-read.table("TCGA.RNAseqV2.lusc.Diff.txt",head=T,sep="\t",as.is=T)
exp.data.sig<-subset(exp.data,Pvalue.ttest<2.4e-06)
head(methy.data)
head(exp.data)

gene.meth<-methy.data[,3]
gene.exp<-sapply(rownames(exp.data.sig),function(x) unlist(strsplit(x,"[|]"))[1])
gene.share<-unique(gene.meth[gene.meth %in% gene.exp])

CHR<-Dmr.Start<-Dmr.End<-Delta.beta<-Meth.Status<-meth.Pvalue<-c()
Gene.Symbol<-Expression.FC<-Expression.Pvalue<-Expression.Status<-c()

for(i in 1:length(gene.share)){
  exp.tmp<-exp.data.sig[which(gene.exp==gene.share[i]),]
  meth.tmp<-methy.data[which(methy.data$gene==gene.share[i]),]
  for(j in 1:nrow(meth.tmp)){
    if(exp.tmp$delta*meth.tmp[j,]$deltaBeta>0){
      CHR<-c(CHR,meth.tmp[j,]$CHR)
      Dmr.Start<-c(Dmr.Start,meth.tmp[j,]$dmr.start)
      Dmr.End<-c(Dmr.End,meth.tmp[j,]$dmr.end)
      Delta.beta<-c(Delta.beta,meth.tmp[j,]$deltaBeta)
      meth.Pvalue<-c(meth.Pvalue,meth.tmp[j,]$dmr.p)
      if(meth.tmp[j,]$deltaBeta>=0){
        Meth.Status<-c(Meth.Status,"Hypomethylation")  
      }else{
        Meth.Status<-c(Meth.Status,"Hypermethylation")  
      }
      Gene.Symbol<-c(Gene.Symbol,gene.share[i])
      fc.tmp<-exp.tmp$mean.A/exp.tmp$mean.B
      Expression.FC<-c(Expression.FC,round(fc.tmp,2))
      Expression.Pvalue<-c(Expression.Pvalue,exp.tmp$Pvalue.ttest)
      if(fc.tmp>=1){
        Expression.Status<-c(Expression.Status,"Overexpression")
      }else{
        Expression.Status<-c(Expression.Status,"lowexpression")
      }
    }
  }
}

output<-data.frame(CHR,Dmr.Start,Dmr.End,Delta.beta,meth.Pvalue,Meth.Status,Gene.Symbol,Expression.FC,Expression.Pvalue,Expression.Status)
dim(output)
head(output)
subset<-subset(output,abs(Delta.beta)>0.2)
write.table(subset,file="TCGA.lusc.Methylation.Expression.Significant.txt",col.names=T,row.names=F,sep="\t",quote=F)
AHRR<-output[grep("AHRR",output$Gene.Symbol),]

#######################################################################################################################
#######################################  Methylation profile plot ##############################################
#######################################################################################################################


load("TCGA.Methy450.lusc.RM.Impute.RData")
library("ggplot2")
map<-methdata$map
data<-methdata$matrix

# situation 1, plot all the DMR genes
DMR<-read.table("TCGA.Meth450.lusc.ChAMP.DMR.txt",head=T,sep="\t",as.is=T)
DMR<-subset(DMR,abs(deltaBeta)>0.1)
dmr.gene<-unique(DMR$gene)

# situation 1, plot some candidate DMR genes
# target.file<-read.table("target.txt",head=T,sep="\t",as.is=T)
# dmr.gene<-unique(target.file$Gene.Symbol)


system("mkdir DMRFigure")
for(i in 1:length(dmr.gene)){  # length(dmr.gene)
  target.gene=dmr.gene[i]
  output.figure.pdf=paste("./DMRFigure/",target.gene,".meth450.pdf",sep="")
  map.tmp<-map[which(target.gene==map[,2]),]
  type=unlist(substr(colnames(data),14,15))
  compare.group=c("01","11")
  x1<-which(type==compare.group[1])
  x2<-which(type==compare.group[2])
  
  if(nrow(map.tmp)>2){
    map.tmp<-map.tmp[order(map.tmp$Chromosome,map.tmp$Genomic_Coordinate),]
    data.tmp<-data[match(map.tmp$Composite.Element.REF,map$Composite.Element.REF),]
    mean1<-apply(data.tmp,1,function(x) mean(x[x1]))
    sd1<-apply(data.tmp,1,function(x) sd(x[x1]))
    mean2<-apply(data.tmp,1,function(x) mean(x[x2]))
    sd2<-apply(data.tmp,1,function(x) sd(x[x2]))
    tmp<-rbind(data.frame(cpg=names(mean1),mean=mean1,sd=sd1,type="Cancer"),data.frame(cpg=names(mean2),mean=mean2,sd=sd2,type="Normal"))
    pdf(output.figure.pdf)
    qplotrlt<-qplot(cpg, mean, colour=type, data=tmp) + 
      xlab("CpG Sites")+
      ylab("Methylation Level (Beta)")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      geom_line(aes(x=cpg,y=mean,group=type))+ 
      geom_errorbar(aes(ymin = mean - 1.96*sd/sqrt(length(x1)), ymax = mean + 1.96*sd/sqrt(length(x2))))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(qplotrlt)
    dev.off()
    print(i)
  }
}

#######################################################################################################################
#######################################  Gene expression plot ##############################################
#######################################################################################################################

expDataSave1="/home/sguo/tcga/lusc/TCGA.RNAseqV2.lusc.RData"    # save RNAseqV2 matrix
exptxtSave1="/home/sguo/tcga/lusc/TCGA.RNAseqV2.lusc.Diff.txt"  # save signficant differential gene ouput
load(expDataSave1)

# situation 1, plot all DGE genes
DGE<-read.table(exptxtSave1,head=T,sep="\t",as.is=T)
DGE<-subset(DGE,Pvalue.ttest<3e-6 & Pvalue.wilcox<3e-6)
dge.gene<-unique(rownames(DGE))

# situation 1, plot some candidate DGE genes
# target.file<-read.table("target.txt",head=T,sep="\t",as.is=T)
# dge.gene<-unique(target.file$Gene.Symbol)


system("mkdir DGEFigure")
library("ggplot2")
for(i in 1:length(dge.gene)){
  target.gene=dge.gene[i]
  type=unlist(substr(colnames(expdata),14,15))
  compare.group=c("01","11")
  x1<-which(type==compare.group[1])
  x2<-which(type==compare.group[2])
  
  match=grep(target.gene,rownames(expdata))
  for(j in 1:length(match)){
    expression<-expdata[match[j],c(x1,x2)]
    output.figure.pdf=paste("./DGEFigure/",gsub("[|]",".",rownames(expdata)[match[j]]),".RNAseqV2.pdf",sep="")
    type<-c(rep("Cancer",length(x1)),rep("Normal",length(x2)))
    data.tmp<-data.frame(expression,type)
    
    pdf(output.figure.pdf)
    qplotrlt<-qplot(type, expression, data=data.tmp, geom=c("boxplot", "jitter"), 
                    fill=type, xlab="", ylab="Gene expression")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    print(qplotrlt)
    dev.off()
    print(i)
  }
}


