
LD<-function(vector){
  rlt<-list()
  table<-matrix(table(vector),2,2)
  pAB=table[1,1]/sum(table)
  pA=(2*table[1,1]+table[2,1]+table[1,2])/(2*sum(table))
  pB=(2*table[2,2]+table[2,1]+table[1,2])/(2*sum(table))
  pa=1-pA
  pb=1-pB
  D=pAB-pA*pB
  if(D>0){
    Dmax=min(pA*pb,pa*pB)
  } else{
    Dmax=max(-pA*pB,-pa*pb)
  } 
  Dp=D/Dmax
  r=Dp/sqrt(pA*pa*pB*pb)
  test<-chisq.test(table)
  chisq<-test$statistic
  
  A1<-as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))
  A2<-as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))  
  fit<-cor.test(A1,A2)
  corr=fit$estimate
  rlt$nobs<-sum(table)
  rlt$corr=as.numeric(corr)
  rlt$corr.p=fit$p.value
  rlt$Dp=Dp
  rlt$r=r
  rlt$p=test$p.value
  return(rlt)
}
