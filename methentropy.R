# methylation entropy
library("entropy")
methentropy<-function(hapinfo){
entropy<-entropy(table(hapinfo),unit="log2")/nchar(hapinfo[1])
return(entropy)
}

