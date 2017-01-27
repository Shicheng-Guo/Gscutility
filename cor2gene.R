cor2gene<-function(cor,refbed){
    bed<-cor2bed(cor)
    rlt<-Rbedtools(functionstring="intersectBed",bed,refbed,opt.string="-wao")
    return(rlt)
}
