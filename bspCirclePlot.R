###   BSP circle plot base on matrix
M<-matrix(sample(c(0,1,1,1),1000,replace=T),10,10)
col=colorRampPalette(c("white", "red"))(20)
circle=c(1,19)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)),ylim=c(0,nrow(M)))
for(i in 1:nrow(M)){
  for(j in 1:ncol(M)){
    points(i,j,col=1,pch=circle[M[i,j]+1],cex=1.5)
    
  }
}


