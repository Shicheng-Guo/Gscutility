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
 rlt$corr=fit$estimate
 rlt$corr.p=fit$p.value
 rlt$p=test$p.value
 rlt$Dp=Dp
 rlt$corr=corr
 rlt$nobs<-sum(table)
 return(rlt)
 }
 vector<-sample(c(rep("CC",3),rep("CT",1),rep("TC",1),rep("TT",3)),100,replace=T)
 LD(vector)
 vector<-sample(c(rep("CC",1),rep("CT",1),rep("TC",1),rep("TT",1)),100,replace=T)
 LD(vector)
 vector<-sample(c(rep("CC",1),rep("CT",2),rep("TC",2),rep("TT",1)),100,replace=T)
 LD(vector)
 mlc<-c()
 for(i in 1:100){
 vector<-sample(c(rep("CC",2),rep("CT",1),rep("TC",1),rep("TT",2)),100,replace=T)
 A1<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))-2))/length(vector)
 A2<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))-2))/length(vector)
 tmp<-data.frame(A1,A2)
 mlc<-rbind(mlc,tmp)
 }
 cor.test(mlc[,1],mlc[,2])
