#!/usr/bin/env Rscript
# For PCA Analysis to methylation 450K dataset
# for ips methylatin 450K analysis
# Oct/4/2016

##########################################################################################################################################
##########################################################################################################################################
library("impute")
library("optparse")
option_list = list(
  make_option(c("-id", "--input"), type="character", default=NULL, help="GSE ID", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

library("GEOquery")
PCAPlot<-function(data,pheno,output,multifigure=T){
  print("PCA Analysis: N rows (objects) x p columns (variables), be sure rownames should be unique")
  pca <- prcomp(data,center=T,scale = F)  
  outputfile=paste(output,".pdf",sep="")
  pdf(outputfile)
  if(multifigure){
    par(mfrow=c(2,2),mar=c(4,4,4,4)) 
  }
  plot((pca$sdev[1:10])^2,type="o",xaxt="n",ylab="Variances",xlab="Principle Components",col="red",lwd=2)
  axis(1,at=0:10,labels=paste("PC",0:10,sep=""))
  var<-c()
  for(i in 1:length(pca$sdev)){var[i]<-sum((pca$sdev[1:i])^2)/sum((pca$sdev)^2)}
  plot(var,ylab="total variance",xlab="number of principle components",lwd=2,type="l")
  abline(h=0.8,col="grey",lty=2)
  abline(v=which(var>0.8)[1],col="grey",lty=2)
  scores <- data.frame(pheno, pca$x[,1:3])
  col = as.numeric(as.factor(pheno))
  plot(x=scores$PC1,y=scores$PC2, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC2),max(scores$PC2)),type="n",xlab="PC1",ylab="PC2")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC2[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.8,lwd=2)
  }
  legend("topleft",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n",pt.lwd=2,,cex=0.5)
  plot(x=scores$PC1,y=scores$PC3, xlim=c(min(scores$PC1),max(scores$PC1)),ylim=c(min(scores$PC3),max(scores$PC3)),type="n",xlab="PC1",ylab="PC3")
  for(i in 1:length(scores$PC1)){
    points(scores$PC1[i],scores$PC3[i],pch=as.numeric(as.factor(pheno))[i],col=col[i],cex=0.9,lwd=2)
  }
  legend("bottomleft",legend=names(table(pheno)),pch=1:length(table(pheno)),col=1:length(table(pheno)),bty="n",pt.lwd=2,cex=0.5)
  dev.off()
}

MDSPlot<-function(mydata,phen,prefix="MDS"){
  print("Classical MDS: N rows (objects) x p columns (variables)[same as PCA], be sure rownames should be unique")
  rownames(mydata)=phen
  output=paste(prefix,".mds.pdf",sep="")
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
ClusterAnalysisPlot<-function(data,phen,prefix="pvclust"){
  colnames(data)=phen
  d <- dist(data,method = "euclidean") # distance matrix
  output=paste(prefix,"euclidean.clust.pdf",sep=".")
  fit <- hclust(d, method="ward") 
  pdf(output)
  plot(fit)
  groups <- cutree(fit, k=6) # cut tree into 5 clusters
  # draw dendogram with red borders around the 5 clusters 
  rect.hclust(fit, k=6, border="red")
  dev.off()
}
PvclusterAnalysisPlot<-function(data,phen,prefix="pvclust"){
  library(pvclust)
  colnames(data)=phen
  fit <- pvclust(data, method.hclust="ward",method.dist="euclidean")
  pdf(paste(prefix,"PvclustAnalysis.pdf",sep="."))
  plot(fit,cex=0.8) # dendogram with p values
  # add rectangles around groups highly supported by the data
  pvrect(fit, alpha=.95)
  dev.off()
}
RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    dat<-data[-NaRAW,]
  }else{
    dat<-data;
  }
  dat
} 

logit<-function (p, percents = range.p[2] > 1, adjust) {
  range.p <- range(p, na.rm = TRUE)
  if (percents) {
    if (range.p[1] < 0 || range.p[1] > 100) 
      stop("p must be in the range 0 to 100")
    p <- p/100
    range.p <- range.p/100
  }
  else if (range.p[1] < 0 || range.p[1] > 1) 
    stop("p must be in the range 0 to 1")
  a <- if (missing(adjust)) {
    if (isTRUE(all.equal(range.p[1], 0)) || isTRUE(all.equal(range.p[2], 
                                                             1))) 
      0.025
    else 0
  }
  else adjust
  if (missing(adjust) && a != 0) 
    warning(paste("proportions remapped to (", a, ", ", 1 - 
                    a, ")", sep = ""))
  a <- 1 - 2 * a
  log((0.5 + a * (p - 0.5))/(1 - (0.5 + a * (p - 0.5))))
}

BMIQ <-function(beta.v,design.v,nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1){
  # Beta MIxture Quantile dilation
  require(RPMM);
  type1.idx <- which(design.v==1);
  type2.idx <- which(design.v==2);
  beta1.v <- beta.v[type1.idx];
  beta2.v <- beta.v[type2.idx];
  ### check if there are exact 0's or 1's. If so, regularise using minimum positive and maximum below 1 values.
  if(min(beta1.v)==0){
    beta1.v[beta1.v==0] <- min(setdiff(beta1.v,0));
  }
  if(min(beta2.v)==0){
    beta2.v[beta2.v==0] <- min(setdiff(beta2.v,0));
  }
  if(max(beta1.v)==1){
    beta1.v[beta1.v==1] <- max(setdiff(beta1.v,1));
  }
  if(max(beta2.v)==1){
    beta2.v[beta2.v==1] <- max(setdiff(beta2.v,1));
  }
  ### estimate initial weight matrix from type1 distribution
  w0.m <- matrix(0,nrow=length(beta1.v),ncol=nL);
  w0.m[which(beta1.v <= th1.v[1]),1] <- 1;
  w0.m[intersect(which(beta1.v > th1.v[1]),which(beta1.v <= th1.v[2])),2] <- 1;
  w0.m[which(beta1.v > th1.v[2]),3] <- 1;
  ### fit type1
  print("Fitting EM beta mixture to type1 probes");
  rand.idx <- sample(1:length(beta1.v),nfit,replace=FALSE)
  em1.o <- blc(matrix(beta1.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
  subsetclass1.v <- apply(em1.o$w,1,which.max);
  subsetth1.v <- c(mean(c(max(beta1.v[rand.idx[subsetclass1.v==1]]),min(beta1.v[rand.idx[subsetclass1.v==2]]))),mean(c(max(beta1.v[rand.idx[subsetclass1.v==2]]),min(beta1.v[rand.idx[subsetclass1.v==3]]))));
  class1.v <- rep(2,length(beta1.v));
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1;
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3;
  nth1.v <- subsetth1.v;
  print("Done");
  ### generate plot from estimated mixture
  if(plots){
    print("Check");
    tmpL.v <- as.vector(rmultinom(1:nL,length(beta1.v),prob=em1.o$eta));
    tmpB.v <- vector();
    for(l in 1:nL){
      tmpB.v <- c(tmpB.v,rbeta(tmpL.v[l],em1.o$a[l,1],em1.o$b[l,1]));
    }
    pdf(paste(sampleID,".Type1fit",".pdf",sep=""),width=6,height=4);
    plot(density(beta1.v));
    d.o <- density(tmpB.v);
    points(d.o$x,d.o$y,col="green",type="l")
    legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
    dev.off();
  }
  ### Estimate Modes 
  d1U.o <- density(beta1.v[class1.v==1])
  d1M.o <- density(beta1.v[class1.v==3])
  mod1U <- d1U.o$x[which.max(d1U.o$y)]
  mod1M <- d1M.o$x[which.max(d1M.o$y)]
  d2U.o <- density(beta2.v[which(beta2.v<0.4)]);
  d2M.o <- density(beta2.v[which(beta2.v>0.6)]);
  mod2U <- d2U.o$x[which.max(d2U.o$y)]
  mod2M <- d2M.o$x[which.max(d2M.o$y)]
  ### now deal with type2 fit
  th2.v <- vector();
  th2.v[1] <- nth1.v[1] + (mod2U-mod1U);
  th2.v[2] <- nth1.v[2] + (mod2M-mod1M);
  ### estimate initial weight matrix 
  w0.m <- matrix(0,nrow=length(beta2.v),ncol=nL);
  w0.m[which(beta2.v <= th2.v[1]),1] <- 1;
  w0.m[intersect(which(beta2.v > th2.v[1]),which(beta2.v <= th2.v[2])),2] <- 1;
  w0.m[which(beta2.v > th2.v[2]),3] <- 1;
  print("Fitting EM beta mixture to type2 probes");
  rand.idx <- sample(1:length(beta1.v),nfit,replace=FALSE)
  em2.o <- blc(matrix(beta2.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
  print("Done");
  ### for type II probes assign to state (unmethylated, hemi or full methylation)
  subsetclass2.v <- apply(em2.o$w,1,which.max);
  subsetth2.v <- c(mean(max(beta2.v[rand.idx[subsetclass2.v==1]]),min(beta2.v[rand.idx[subsetclass2.v==2]])),mean(max(beta2.v[rand.idx[subsetclass2.v==2]]),min(beta2.v[rand.idx[subsetclass2.v==3]])));
  class2.v <- rep(2,length(beta2.v));
  class2.v[which(beta2.v < subsetth2.v[1])] <- 1;
  class2.v[which(beta2.v > subsetth2.v[2])] <- 3;
  ### generate plot
  if(plots){
    tmpL.v <- as.vector(rmultinom(1:nL,length(beta2.v),prob=em2.o$eta));
    tmpB.v <- vector();
    for(lt in 1:nL){
      tmpB.v <- c(tmpB.v,rbeta(tmpL.v[lt],em2.o$a[lt,1],em2.o$b[lt,1]));
    }
    pdf(paste(sampleID,".Type2fit",".pdf",sep=""),width=6,height=4);
    plot(density(beta2.v));
    d.o <- density(tmpB.v);
    points(d.o$x,d.o$y,col="green",type="l")
    legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
    dev.off();
  }
  classAV1.v <- vector();classAV2.v <- vector();
  for(l in 1:nL){
    classAV1.v[l] <-  em1.o$mu[l,1];
    classAV2.v[l] <-  em2.o$mu[l,1];
  }
  
  ### start normalising type2 probes
  print("Start normalising type 2 probes");
  nbeta2.v <- beta2.v;
  ### select U probes
  lt <- 1;
  selU.idx <- which(class2.v==lt);
  selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])];
  selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])];
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selUR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
  nbeta2.v[selUR.idx] <- q.v;
  p.v <- pbeta(beta2.v[selUL.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=TRUE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=TRUE);
  nbeta2.v[selUL.idx] <- q.v;
  ### select M probes
  lt <- 3;
  selM.idx <- which(class2.v==lt);
  selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])];
  selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])];
  ### find prob according to typeII distribution
  p.v <- pbeta(beta2.v[selMR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
  ### find corresponding quantile in type I distribution
  q.v <- qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
  nbeta2.v[selMR.idx] <- q.v;
  if(doH){ ### if TRUE also correct type2 hemimethylated probes
    ### select H probes and include ML probes (left ML tail is not well described by a beta-distribution).
    lt <- 2;
    selH.idx <- c(which(class2.v==lt),selML.idx);
    minH <- min(beta2.v[selH.idx])
    maxH <- max(beta2.v[selH.idx])
    deltaH <- maxH - minH;
    #### need to do some patching
    deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
    deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])
    ## new maximum of H probes should be
    nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM;
    ## new minimum of H probes should be
    nminH <- max(nbeta2.v[selU.idx]) + deltaUH;
    ndeltaH <- nmaxH - nminH;
    ### perform conformal transformation (shift+dilation)
    ## new_beta_H(i) = a + hf*(beta_H(i)-minH);
    hf <- ndeltaH/deltaH ;
    ### fix lower point first
    nbeta2.v[selH.idx] <- nminH + hf*(beta2.v[selH.idx]-minH);
  }
  pnbeta.v <- beta.v;
  pnbeta.v[type1.idx] <- beta1.v;
  pnbeta.v[type2.idx] <- nbeta2.v;
  ### generate final plot to check normalisation
  if(plots){
    print("Generatingfinal plot");
    d1.o <- density(beta1.v);
    d2.o <- density(beta2.v);
    d2n.o <- density(nbeta2.v);
    ymax <- max(d2.o$y,d1.o$y,d2n.o$y);
    pdf(paste(sampleID,"CheckBMIQ",".pdf",sep=""),width=6,height=4)
    plot(density(beta2.v),type="l",ylim=c(0,ymax),xlim=c(0,1));
    points(d1.o$x,d1.o$y,col="red",type="l");
    points(d2n.o$x,d2n.o$y,col="blue",type="l");
    legend(x=0.5,y=ymax,legend=c("type1","type2","type2-BMIQ"),bty="n",fill=c("red","black","blue"));
    dev.off();
  }
  print(paste("Finished for sample ",sampleID,sep=""));
  return(list(nbeta=pnbeta.v,class1=class1.v,class2=class2.v,av1=classAV1.v,av2=classAV2.v,hf=hf,th1=nth1.v,th2=th2.v));
}
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
  chart = ggplot(data.frame(tsne_iris), aes(xtsne,ytsne))+geom_point(size=1,alpha=1,aes(colour = factor(phen1)))+ggtitle("tSNE dimensions colored by digit")
  # change the size and alpha could make great beautful figure
  chart
  dev.off()
  
}

##########################################################################################################################################
##########################################################################################################################################

GSEID<-opt$input  # GSEID="GSE51954"
GEO <- getGEO(GSEID,destdir=getwd())
library("GEOquery")
save(GEO,file=paste(GSEID,".RData",sep=""))
load(paste(GSEID,"_matrix.Rdata",sep=""))
beta <- as.data.frame(exprs(GEO[[1]]))
phen <- pData(phenoData(GEO[[1]]))

#########################
phen1<-unlist(lapply(as.character(phen$characteristics_ch1),function(x) unlist(strsplit(x,"[: ]"))[length(unlist(strsplit(x,"[: ]")))]))
phen2<-unlist(lapply(as.character(phen$characteristics_ch1.1),function(x) unlist(strsplit(x,"[: ]"))[length(unlist(strsplit(x,"[: ]")))]))
phen3<-unlist(lapply(as.character(phen$characteristics_ch1.2),function(x) unlist(strsplit(x,"[: ]"))[length(unlist(strsplit(x,"[: ]")))]))
phen4<-unlist(lapply(as.character(phen$characteristics_ch1.3),function(x) unlist(strsplit(x,"[: ]"))[length(unlist(strsplit(x,"[: ]")))]))

table(phen1)
table(phen2)
table(phen3)
#########################

# modify
# phen3<-as.numeric(unlist(lapply(phen3,function(x) gsub("y","",x))))

################################################################
########## assess the imputation accuracy #######################
#################################################################
library("impute")
BETA<-as.numeric(data.matrix(beta))
Sample<-sample(1:length(BETA),20000)
Value<-BETA[Sample]
BETA[Sample]<-NA
BETA2<-matrix(BETA,ncol=ncol(beta),byrow=F)
BETA.normal<-impute.knn(data.matrix(BETA2),rowmax = 0.5)
BETA.normal2<-as.numeric(BETA.normal$data)
Value2<-BETA.normal2[Sample]
Error<-Value-Value2

fileoutput<-paste(GSEID,".knn.imputation.error.pdf",sep="")
pdf(fileoutput)
plot(hist(Error,col="blue",xlab="Real-impute",ylab="Number"))
dev.off()

#############################################################################
########## assess the quantile normalization accuracy #######################
#############################################################################
library("impute")
library("preprocessCore")
beta.normal<-impute.knn(data.matrix(beta),rowmax = 0.8)
Beta=na.omit(beta.normal$data)
dim(Beta)
# NewBeta<-normalize.quantiles(Beta,copy=TRUE)
# rownames(NewBeta)=rownames(Beta)
# fileoutput<-paste(GSEID,".quantile.error.pdf",sep="")
# pdf(fileoutput)
# plot(hist(Beta-NewBeta,col="blue",xlab="Real-impute",ylab="Number"))
# dev.off()

mapinfo<-read.table("/home/shg047/oasis/db/hg19/GPL13534.sort.bed")
design<-as.numeric(mapinfo[match(rownames(Beta),mapinfo[,4]),7])
NewBetaBMIQ<-BMIQ(Beta,design,sampleID=GSEID)
fileoutput<-paste(GSEID,".quantile.bmiq.error.pdf",sep="")
pdf(fileoutput)
plot(hist(Beta-NewBetaBMIQ$nbeta,col="blue",xlab="Real-impute",ylab="Number"))
dev.off()

############################################################################
################## PCA Analysis ########################################
############################################################################
fileoutput<-paste(GSEID,"pca.phen1.pdf",sep=".")
PCAPlot(t(NewBetaBMIQ$nbeta),phen1,output=fileoutput,multifigure=T)
fileoutput<-paste(GSEID,"pca.phen2.pdf",sep=".")
PCAPlot(t(NewBetaBMIQ$nbeta),phen2,output=fileoutput,multifigure=T)
fileoutput<-paste(GSEID,"pca.phen3.pdf",sep=".")
PCAPlot(t(NewBetaBMIQ$nbeta),phen3,output=fileoutput,multifigure=T)
fileoutput<-paste(GSEID,"pca.phen4.pdf",sep=".")
PCAPlot(t(NewBetaBMIQ$nbeta),phen4,output=fileoutput,multifigure=T)
############################################################################
###################### TSNE Analysis #######################################
############################################################################
library("ggplot2")
HVFinput<-NewBetaBMIQ$nbeta[order(unlist(apply(NewBetaBMIQ$nbeta,1,sd)),decreasing=T)[1:2000],]

fileoutput<-paste(GSEID,prefix="tsne.phen1",sep=".")
TSNEAnalysis(t(HVFinput),phen1,prefix=fileoutput)

fileoutput<-paste(GSEID,prefix="tsne.phen2",sep=".")
TSNEAnalysis(t(HVFinput),phen2,prefix=fileoutput)

fileoutput<-paste(GSEID,prefix="tsne.phen3",sep=".")
TSNEAnalysis(t(HVFinput),phen3,prefix=fileoutput)

fileoutput<-paste(GSEID,prefix="tsne.phen4",sep=".")
TSNEAnalysis(t(HVFinput),phen4,prefix=fileoutput)

############################################################################
################## Cluster Analysis ########################################
############################################################################
print("classic cluster analysis.....");
outputfile=paste(GSEID,"cluster.pdf",sep=".")
ClusterAnalysisPlot(HVFinput,phen1,prefix=outputfile)
outputfile=paste(GSEID,"cluster.pdf",sep=".")
ClusterAnalysisPlot(HVFinput,phen2,prefix=outputfile)
outputfile=paste(GSEID,"cluster.pdf",sep=".")
ClusterAnalysisPlot(HVFinput,phen3,prefix=outputfile)
outputfile=paste(GSEID,"cluster.pdf",sep=".")
ClusterAnalysisPlot(HVFinput,phen4,prefix=outputfile)

print("pvcluster analysis.....");
outputfile=paste(GSEID,"pvcluster.phen1.pdf",sep=".")
PvclusterAnalysisPlot(HVFinput,phen1,prefix=outputfile)
outputfile=paste(GSEID,"pvcluster.phen2.pdf",sep=".")
PvclusterAnalysisPlot(HVFinput,phen2,prefix=outputfile)
outputfile=paste(GSEID,"pvcluster.phen3.pdf",sep=".")
PvclusterAnalysisPlot(HVFinput,phen3,prefix=outputfile)
outputfile=paste(GSEID,"pvcluster.phen4.pdf",sep=".")
PvclusterAnalysisPlot(HVFinput,phen4,prefix=outputfile)

############################################################################
################## MDS analysis ########################################
############################################################################
print("MDS analysis.....");
prefix=paste(GSEID,"mds.phen1",sep=".")
MDSPlot(t(NewBetaBMIQ$nbeta),phen=phen1,prefix=prefix)
prefix=paste(GSEID,"mds.phen2",sep=".")
MDSPlot(t(NewBetaBMIQ$nbeta),phen=phen2,prefix=prefix)
prefix=paste(GSEID,"mds.phen3",sep=".")
MDSPlot(t(NewBetaBMIQ$nbeta),phen=phen3,prefix=prefix)
prefix=paste(GSEID,"mds.phen4",sep=".")
MDSPlot(t(NewBetaBMIQ$nbeta),phen=phen4,prefix=prefix)

############################################################################
################## Differential Analysis ###################################
############################################################################


############################################################################
######################## CCA Analysis ######################################
############################################################################

