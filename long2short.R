long2short<-function(barcode){
  return(as.array(str_extract(barcode,"TCGA-[0-9|a-z|A-Z]*-[0-9|a-z|A-Z]*-[0-9]*")))
}
