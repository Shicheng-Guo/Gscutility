setwd("/home/sguo/mice/scnt12vips")
library("DSS")
read_methBED <- function(filename)
{
   meth = read.table(file=filename,stringsAsFactors=F,skip=1)
   x = meth[,c(1,2,5)]
   x = cbind(x,as.integer(meth[,5]*meth[,4]))
   colnames(x) = c('chr','pos','N','X');
   return(x)
}
############################################################
#Parameters
fdr_cutoff = 0.05

#Input
filenames_1 = c(
            'miPS_B3.BED.txt.trim',
            'miPS_1E12P20.BED.txt.trim','miPS_2A4F1.BED.txt.trim','miPS_2A4F33.BED.txt.trim'
            )
samplenames_1 = c('miPS2','miPS1E12','miPS2A1','miPS2A33')
filenames_2 = c(
            'SCNT_B12.BED.txt.trim',
            'SCNT_NB3.BED.txt.trim'
            )
samplenames_2 = c('SCNT_B12','SCNT_NB3')
############################################################
#Read input files
input = list()
filenames=c(filenames_1,filenames_2)
for(file in filenames)
{
  input = append(input,list(read_methBED(file)))
}
save(input,file="input.RData")
#Calling differentially methylated sites
BSobj <- makeBSseqData(input,c("miPS2","miPS1E12","miPS2A1","miPS2A33", "SCNT_B12", "SCNT_NB3") )
save(BSobj,file="BSobj.RData")
dmlTest <- DMLtest(BSobj, group1=c("miPS2", "miPS1E12","miPS2A1","miPS2A33"), group2=c("SCNT_B12", "SCNT_NB3"),smoothing = T,smoothing.method = c("BSmooth"), smoothing.span = 200)
save(dmlTest,file="dmlTest.region.RData")
DMS = callDML(dmlTest)
#Write output to hard disk
fdr_cutoff=0.05
index = which(DMS$fdr <= fdr_cutoff)
write.table(DMS[index,],file="output_DMS_4iPS_SCNTB12NB3.region.tsv",row.names=F,quote=F,sep="\t")
