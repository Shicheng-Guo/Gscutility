chan<-function(x){
  x[x<0.5]<-0; x[x>=0.5]<-1
  x
}
assess<-function(model,train,test,trainlabel,testlabel){
  t1<-table(chan(predict(model, train)), trainlabel)
  t2<-table(chan(predict(model, test)), testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}
assess2<-function(model1,model2,trainlabel,testlabel){
  t1<-table(model1,trainlabel)
  t2<-table(model2,testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}
assess3<-function(model,train,test,trainlabel,testlabel){
  t1<-table((predict(model, train)), trainlabel)
  t2<-table((predict(model, test)), testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}
assess4<-function(model1,model2,trainlabel,testlabel){
  t1<-table(chan(model1), trainlabel)
  t2<-table(chan(model2), testlabel)
  sen.train<-t1[2,2]/(t1[1,2]+t1[2,2])
  spe.train<-t1[1,1]/(t1[2,1]+t1[1,1])
  sen.test<-t2[2,2]/(t2[1,2]+t2[2,2])
  spe.test<-t2[1,1]/(t2[2,1]+t2[1,1])
  accu.train<-(t1[1,1]+t1[2,2])/sum(t1)
  accu.test<-(t2[1,1]+t2[2,2])/sum(t2)
  rlt<-c(sen.train,spe.train, accu.train,sen.test,spe.test,accu.test)
  rlt
}
