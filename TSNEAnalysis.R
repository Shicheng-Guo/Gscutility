TSNEAnalysis<-function(mydata,phen,prefix="TSNE"){
  print("t-sne analysis: N rows (objects) x P columns (variables) [same as PCA], be sure rownames should be unique")
  library("tsne")
  output=paste(prefix,".tsne-1.pdf",sep="")
  pdf(output)
  colors = rainbow(length(unique(phen)))
  names(colors) = unique(phen)
  ecb = function(x,y){ plot(x,t='n',xlab="Coordinate 1",ylab="Coordinate 2"); text(x,labels=phen, col=colors[phen]) }
  tsne_iris = data.frame(tsne(mydata, epoch_callback = ecb, perplexity=50))
  dev.off()
  
  colnames(tsne_iris)<-c("xtsne","ytsne")
  output=paste(prefix,".tsne-2.pdf",sep="")
  pdf(output)
  chart = ggplot(data.frame(tsne_rlt), aes(xtsne,ytsne))+geom_point(size=1,alpha=1,aes(colour = factor(phen1)))+ggtitle("tSNE dimensions colored by digit")
  # change the size and alpha could make great beautful figure
  chart
  dev.off()
  
}
