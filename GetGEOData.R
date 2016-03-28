#!/usr/bin/env Rscript
library("optparse")
Dir=getwd()
option_list = list(
  make_option(c("-i", "--GSEID"), type="character", default=NULL, help="dataset file name", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$GSEID)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

library("GEOquery")
GEOSet <- getGEO(opt$GSEID)
data <- as.data.frame(exprs(GEOSet[[1]]))
phen <- pData(phenoData(GEOSet[[1]]))
fileName1<-paste(opt$GSEID,"_matrix.Rdata",sep="")
fileName2<-paste(opt$GSEID,"_matrix.phen.txt",sep="")
fileName3<-paste(opt$GSEID,"_matrix.data.txt",sep="")
GEO<-list()
rm(GEOSet)
GEO$data<-data
GEO$phen<-phen
save(GEO, file=fileName1)
write.table(GEO$phen,file=fileName2,sep="\t",quote=F)
write.table(GEO$data,file=fileName3,sep="\t",quote=F)
