###   Figure LD plot, heatmap for LD blocks
library("grDevices")
col=colorRampPalette(c("white", "red"))(20)
M <- matrix(runif(400),20,20)
M[lower.tri(M)] <- NA
M[M>0.4]<-1
image(M,col = col,frame=F,xaxt="n",yaxt="n")
