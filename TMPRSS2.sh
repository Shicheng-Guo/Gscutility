####################################################################################################################################
mkdir TMPRSS2
awk '{print $1,$2}' OFS="\t" missense.txt > TMPRSS2.hg19.txt
bcftools view chr21.vcf.gz -R ./TMPRSS2/TMPRSS2.hg19.txt -Oz -o ./TMPRSS2/TMPRSS2.vcf.gz
plink --vcf TMPRSS2.vcf.gz --make-bed --out TMPRSS2

cp ../../phase3_corrected.psam ./
awk '{print $1,$1,$4,$5,$6}' phase3_corrected.psam > phase3.psam

for i in ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI 
do
grep $i phase3.psam > $i.txt
plink --bfile TMPRSS2 --keep $i.txt --freq --out $i
done

file=list.files(pattern="*.frq")
FRQ<-c()
for(i in 1:length(file)){
  frq<-read.table(file[i],head=T)
  FRQ<-cbind(FRQ,frq[,5])
}
colnames(FRQ)<-unlist(lapply(strsplit(file,"[.frq]"),function(x) x[1]))
rownames(FRQ)<-frq[,2]
write.csv(FRQ,file="TMPRSS2.frq.csv",row.names=T,quote=F)
####################################################################################################################################

