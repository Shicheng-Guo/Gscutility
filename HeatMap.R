
HeatMap<-function(data,phen,plot="heatmap.pdf",cexRow = 0.2,cexCol = 0.9,Colv=T,Rowv=T){
#  subset<-order(unlist(apply(data,1,function(x) sd(x,na.rm=T))),decreasing=T)[1:4000]
  x1<-which(phen==unique(phen)[1])
  x2<-which(phen==unique(phen)[2])
  
  subset<-order(unlist(apply(data,1,function(x) t.test(x[x1],x[x2])$p.value)),decreasing=F)[1:4000]  
  data<-data[subset,]
  library("gplots")
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  colors <-greenred(75)
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  ColSideColors=sidecol(phen)
  pdf(plot)
  heatmap.2(data,trace="none",cexRow = cexRow,cexCol = cexCol, ColSideColors=ColSideColors,density.info="none",col=colors,Colv=Colv,Rowv=Rowv,keysize=0.9, margins = c(5, 10))
  dev.off()
}

HeatMap(data=myNorm$beta,phen=c(rep('T',3),rep('C',2),rep('T',7)),plot="heatmap.pdf")


