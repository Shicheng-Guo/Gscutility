write.bed<-function(bed,file){
  if(ncol(bed)==3){
  bed[,4]<-paste(bed[,1],":",bed[,2],"-",bed[,3],sep="")  
  }
  if(ncol(bed)>=4){
    write.table(bed,file=file,sep="\t",col.names=F,row.names=F,quote=F)
  }
}
