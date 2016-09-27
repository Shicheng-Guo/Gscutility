MDSPlot<-function(mydata,phen,prefix="MDS"){
rownames(mydata)=phen
output=paste(prefix,".mds.pdf",sep="")
print("Classical MDS: N rows (objects) x p columns (variables)[same as PCA], be sure rownames should be unique")
d <- dist(mydata) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
pdf(output)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric	MDS",	type="n")
text(x, y, labels = row.names(mydata), cex=.7)
dev.off()
}
