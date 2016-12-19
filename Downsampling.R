#!/usr/bin/R
# Downsampling data by the index (Group ID, Cell type) at defined propotion
# Data: 12/19/2016

# raw data simulation
data<-abs(matrix(rnorm(450),150,3))
idx<-sample(1:4,150,replace=T)
input<-data.frame(data,idx)
dim(input)
head(input)
prop=0.5

# Function for downsampling
downsamplingbyclust<-function(data,idx,prop){
input<-data.frame(data,idx)
Row<-c()
for(i in unique(input$idx)){
 rowc<-sample(which(input$idx==i),round(prop*length(which(input$idx==i)))) 
 Row<-c(Row,rowc)
}
return(Row)
}

# output the downsampling data
down<-input[downsamplingbyclust(data,idx,0.5),]
write.table(down,file="dowsample.txt",row.names=F,col.names=T,sep="\t",quote=F)
