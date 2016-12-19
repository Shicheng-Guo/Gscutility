#!/usr/bin/R
# Downsampling data by the index (Group ID, Cell type) at defined propotion
# for Dr. Blue Lake
# Data: 12/19/2016

# raw data simulation
idx<-sample(1:2,10,replace=T)
idx

# Function for downsampling
downsamplingbyclust<-function(idx,number){
INDEX<-c()
for(i in unique(idx)){
 if(length(which(idx==i))>=number){
   index<-sample(which(idx==i),number) 
 }else{
   index<-which(idx==i)
 }
  INDEX<-c(INDEX,index)
}
return(INDEX)
}

# output the downsampling data
position<-downsamplingbyclust(idx,number=5)
DATA<-data[,position]
write.table(DATA,file="dowsample.txt",row.names=F,col.names=T,sep="\t",quote=F)
