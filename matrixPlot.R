#!/usr/bin/env Rscript
library("optparse")
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="dataset file name", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="out.txt",help="output file name [default= %default]", metavar="character"),
  make_option(c("-xlab", "--xlab"), type="character", default="",help="xlab of the barplot", metavar="character"),
  make_option(c("-ylab", "--ylab"), type="character", default="",help="ylab of the barplot", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}


data<-read.table(opt$input,sep="\t",head=T,row.names=1,check.names=F,as.is=T)
miss.ratio<-sum(is.na(data))/(nrow(data)*ncol(data))
print(paste("Miss ratio of the dataset is:",sprintf("%1.2f%%",100*miss.ratio),sep=" "))

newdata<-data.frame(value=as.numeric(data.matrix(data)),Group=rep(colnames(data),each=nrow(data)))

myData <- aggregate(newdata$value,by =list(type=newdata$Group),
                    FUN = function(x) c(mean = mean(x,na.rm=T),
                                        sd = sd(x,na.rm=T),
                                        sem=sd(x,na.rm=T)/sqrt(length(na.omit(x))),
                                        min=min(x,na.rm=T),
                                        max=max(x,na.rm=T),
                                        median=median(x,na.rm=T),
                                        me=qt(1-0.05/2,df=length(na.omit(x))*sd(x,na.rm=T)/sqrt(length(na.omit(x)))))
)
myData <- do.call(data.frame, myData)
colnames(myData)=c("type","mean","sd","sem","min","max","median","me")
Tmp<-which(myData$sd-myData$mean<0)
myData$sd[Tmp]<-myData$mean[Tmp]

library("ggplot2")
output1<-paste(opt$output,"barplot.pdf",sep=".")
output2<-paste(opt$output,"boxplot.pdf",sep=".")

pdf(output1)
ylab=opt$ylab
xlab=opt$xlab
ggplot(myData, aes(x =type, y = mean)) +
  geom_bar(position = position_dodge(), stat="identity", fill="blue",width=0.8) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),size=1.0) +
  theme_bw() +
  theme(panel.grid.major = element_blank())+
  coord_flip()+
  xlab(xlab) +
  ylab(ylab)+
  theme(axis.text=element_text(),axis.title=element_text(),axis.text.y = element_text(hjust=0))
dev.off()

pdf(output2)
ggplot(newdata, aes(factor(Group),value)) +
  geom_boxplot(aes(fill = factor(Group)),outlier.shape=NA)+
  coord_flip()+
  theme(axis.text=element_text(),axis.title=element_text(),axis.text.y = element_text(hjust=0),legend.position="none")
dev.off()

