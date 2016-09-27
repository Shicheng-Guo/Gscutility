HeatMap<-function(data,phen,varselect=1000,plot="heatmap.pdf",cexRow = 0.01,cexCol = 1.2,Colv=T,Rowv=T){
  print("Heatmap plot: p rows (variables) x N columns (objects), be sure rownames should be unique")
  library("gplots")
  #  subset<-order(unlist(apply(data,1,function(x) sd(x,na.rm=T))),decreasing=T)[1:4000]
  x1<-which(phen==unique(phen)[1])
  x2<-which(phen==unique(phen)[2])
  subset<-order(unlist(apply(data,1,function(x) t.test(x[x1],x[x2])$p.value)),decreasing=F)[1:varselect]  
  data<-data[subset,]
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
