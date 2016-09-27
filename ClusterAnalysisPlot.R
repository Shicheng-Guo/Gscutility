ClusterAnalysisPlot<-function(data,pheno,suffix="pvclust"){
library(pvclust)
fit <- pvclust(data, method.hclust="ward",method.dist="euclidean")
pdf(paste(suffix,"clustAnalysis.pdf",sep="."))
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)
dev.off()
}
