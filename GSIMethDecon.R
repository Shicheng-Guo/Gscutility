GSIMethDecon<-function(data){
  data<-data.matrix(data)
  group=names(table(colnames(data)))
  index=colnames(data)
  gsi<-gmaxgroup<-avebase<-c()
  for(i in 1:nrow(data)){
  gsit<-0
  gmax<-names(which.max(tapply(as.numeric(data[i,]),index,function(x) mean(x,na.rm=T))))
  for(j in 1:length(group)){
    tmp<-(1-10^(mean(data[i,][which(index==group[j])],na.rm=T))/10^(mean(data[i,][which(index==gmax)],,na.rm=T)))/(length(group)-1)
    gsit<-gsit+tmp
  }
  ave<-tapply(data[i,], index, function(x) mean(x,na.rm=T))
  gmaxgroup<-c(gmaxgroup,gmax)
  gsi<-c(gsi,gsit)
  avebase<-rbind(avebase,ave)
  print(i)
  }
  rlt=data.frame(region=rownames(data),group=gmaxgroup,GSI=gsi,AVE=avebase)
  rlt<-gsi[-which(gsi[,2]!="WBC" & gsi[,grep("WBC",colnames(gsi))]>0.1),]
  return(rlt)
}
