bed3tobed4<-function(bed,extend=0){
  bed[,2]<-as.numeric(as.character(bed[,2]))-extend
  bed[,3]<-as.numeric(as.character(bed[,3]))+extend
  bed[,4]<-paste(bed[,1],":",bed[,2],"-",bed[,3],sep="")  
  return(bed)
}
