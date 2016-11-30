###   BSP circle plot base on matrix
M<-matrix(sample(c(sample(0,100,replace=T),1),400,replace=T),100,4)
col=colorRampPalette(c("white", "red"))(20)
circle=c(1,19)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)+1),ylim=c(0,nrow(M)+1))
for(i in 1:ncol(M)){
  for(j in 1:nrow(M)){
    points(i,j,col=1,pch=circle[M[j,i]+1],cex=1)
    
  }
}



# example:

M<-matrix(c(rep(c(0,1,0,1),each=6),c(rep(c(1,0,1,0),each=6))),12,4)
M<-matrix(c(rep(c(0,1,0,1),each=6),c(rep(c(0,1,0,1),each=6))),12,4)
M<-matrix(c(rep(c(0,1),time=6),rep(c(1,0),time=6),rep(c(0,1),time=6),rep(c(1,0),time=6)),12,4)

col=colorRampPalette(c("white", "black"))(20)
circle=c(1,19)
plot(x=nrow(M),y=ncol(M),type="n",xlab="",ylab="",xlim=c(0,ncol(M)+1),ylim=c(0,nrow(M)+1))
for(i in 1:ncol(M)){
  for(j in 1:nrow(M)){
    points(i,j,col=1,pch=circle[M[j,i]+1],cex=1)
    
  }
}
