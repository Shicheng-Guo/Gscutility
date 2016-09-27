TSNEAnalysis<-function(mydata,phen,prefix="TSNE"){
  print("t-sne analysis: N rows (objects) x P columns (variables) [same as PCA], be sure rownames should be unique")
  library("tsne")
  output=paste(prefix,".tsne.pdf",sep="")
  pdf(output)
  colors = rainbow(length(unique(phen)))
  names(colors) = unique(phen)
  ecb = function(x,y){ plot(x,t='n',xlab="Coordinate 1",ylab="Coordinate 2"); text(x,labels=phen, col=colors[phen]) }
  tsne_iris = tsne(mydata, epoch_callback = ecb, perplexity=50)
  dev.off()
}
