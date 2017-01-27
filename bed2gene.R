bed2gene<-function(bed,refbed){
    rlt<-Rbedtools(functionstring="intersectBed",bed,refbed,opt.string="-wao")
    return(rlt)
}
