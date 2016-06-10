
tsneplot<-function(mydata,phen,plot="tsne.plot.pdf"){
library("tsne")
data=data.frame(phen,mydata)
## Not run:
pdf(plot)
colors = rainbow(length(unique(data$phen)))
names(colors) = unique(data$phen)
ecb = function(x,y){plot(x,t='n'); text(x,labels=data$phen, col=colors[data$phen]) }
tsne_iris = tsne(data, epoch_callback = ecb, perplexity=10)
dev.off()
}
tsneplot(mydata[,1:ncol(mydata)],phen,plot="esca.tsne.plot.pdf")
