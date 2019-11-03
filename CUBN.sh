# 2019/11/02
cd /gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R3
input="CUBN"
mkdir $input
grep -w "CUBN" ~/hpc/db/hg19/refGene.hg19.V2.bed | awk '{print $1,$2-5000000,$3+5000000,$5}' OFS="\t" | sort -u | bedtools sort -i > $input.hg19.bed
for i in {1..22} X 
do 
bcftools view -R $input.hg19.bed chr$i.dose.dbSNP.hg19.vcf.gz -Oz -o ./$input/chr$i.dose.dbSNP.miRNA.hg19.vcf.gz
tabix -p vcf ./$input/chr$i.dose.dbSNP.miRNA.hg19.vcf.gz
echo $i
done

cd $input
wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/manhattan.qqplot.R -O manhattan.plot.R
wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/make.fancy.locus.plot.unix.R -O make.fancy.locus.plot.unix.R
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/localhit.pl -O localhit.pl
ls chr*.dose.dbSNP.miRNA.hg19.vcf.gz > filelist.txt
bcftools concat -f filelist.txt -Oz -o ROI.hg19.vcf.gz
plink --vcf ROI.hg19.vcf.gz --make-bed --out ROI.dbsnp
perl ../plink/phen.pl ROI.dbsnp.fam > mylist.txt
plink --vcf ROI.hg19.vcf.gz --keep mylist.txt --make-bed --out ROI.dbsnp
perl ../plink/phen.pl ROI.dbsnp.fam > ROI.dbsnp.new
wc -l ROI.dbsnp.fam
wc -l ROI.dbsnp.new
mv ROI.dbsnp.new ROI.dbsnp.fam
plink --bfile ROI.dbsnp --freq --logistic --adjust --out ROI.dbsnp
plink --bfile ROI.dbsnp --assoc --adjust --out ROI.dbsnp
plink --bfile ROI.dbsnp --counts --assoc --adjust --out ROI.dbsnp.counts

touch readme.txt
bcftools view -G ROI.hg19.vcf.gz -Ov -o ROI.refGene.vcf.txt
echo  "Title: Novel Genetic susceptibility loci in" $input "associated with seropositive rheumatoid arthritis" > readme.txt
echo "" >> readme.txt
plink --bfile ROI.dbsnp --hardy --out ROI
echo  "ROI.hwe: Hardy-Weinberg test statistics for each SNP for" $input. "SNPs in control group should not signficant which means P should higher than 0.05 or 0.01. this table will be put to supplementary table" >> readme.txt
plink --bfile ROI.dbsnp --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI
echo  "ROI.assoc.logistic: logistic based case-control test for" $input. "default style is to test additive model --ADD-- in logistic regression. this file will be one of most important table in the manuscript" >> readme.txt
echo  "ROI.assoc.logistic.adjusted: this file include all the multiple-test corrected P-value for each SNPs in" $input. " When you prepare the manuscript, this file should be integrate with above file --ROI.assoc.logistic--" >> readme.txt
plink --bfile ROI.dbsnp --assoc --counts --adjust --ci 0.95 --out ROI
echo  "ROI.assoc: Chi-square based case-control test for" $input. "this file will be one of most important table in the manuscript since it showed the number of alleles in case and control" >> readme.txt
echo  "ROI.assoc.adjusted: this file include all the multiple-test corrected P-value --in the file: ROI.assoc-- for each SNPs in" $input. " When you prepare the manuscript, this file should be integrate with above file --ROI.assoc--" >> readme.txt
plink --bfile ROI.dbsnp --fisher --counts --adjust --ci 0.95 --out ROI
echo  'ROI.assoc.fisher: Fisher exact test based case-control association between SNPs and RA. This file will be useful when any cell lt 5. Usually when certain cell have number lt 5, we report fisher P-value not Chi-square P-value' >> readme.txt
echo  "ROI.assoc.fisher.adjusted: this file include all the multiple-test corrected P-value --in the file: ROI.assoc-- for each SNPs in" $input. " When you prepare the manuscript, this file should be integrate with above file --ROI.assoc--" >> readme.txt
plink --bfile ROI.dbsnp --model fisher --ci 0.95 --out ROI
echo  "ROI.model: Fisher's exact test based case-control association with different models for" $input. "this file is one of most important table in the manuscript" >> readme.txt
plink --bfile ROI.dbsnp --logistic --dominant --ci 0.95 --out ROI.dominant
echo  "ROI.dominant.assoc.logistic: logistic regression based association in dominant models for" $input. "this file is one of most important table in the manuscript. In the file of ROI.model, the DOM is based on fisher's exact test" >> readme.txt
plink --bfile ROI.dbsnp --logistic --recessive --ci 0.95 --out ROI.recessive
echo  "ROI.recessive.assoc.logistic: logistic regression based association in recessive models for" $input. "this file is one of most important table in the manuscript" >> readme.txt
~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out haplotype --noweb
echo  'haplotype.assoc.hap: chi-square test based haplotype association. This file is important which can be shown with significant haplotype as a table in the manuscript' >> readme.txt
awk '$6=="NMISS"{print}' ROI.assoc.logistic > ROI.assoc.logistic.add
awk '$5=="ADD"{print}' ROI.assoc.logistic >> ROI.assoc.logistic.add

sort -k12n,12 ROI.assoc.logistic | head -n 3
sort -k12n,12 ROI.assoc.logistic | tail -n 3

rs="rs2356827"
chr=10
grep $rs ROI.assoc.logistic
cp ROI.dbsnp.fam $rs.fam
plink --bfile ROI.dbsnp --snp $rs --recode --out $rs
plink --bfile ROI.dbsnp --r2 --ld-snp $rs --ld-window-kb 100 --ld-window 99999 --ld-window-r2 0 --out $rs
perl localhit.pl $rs > $rs.local
Rscript make.fancy.locus.plot.unix.R $rs $rs $chr $rs.local 6 0.05

rm *.R
rm *.pl
echo "" >> readme.txt
echo  'Abstract: The heritability of RA has been shown from twin studies to be 60%. Since 2007, rapid advances in technology underpinning the use of genome-wide association studies have allowed the identification of hundreds of genetic risk factors for many complex diseases. There are now >100 genetic loci that have been associated with RA. In the previous study, the contribution of HLA to heritability has been estimated to be 11–37% while 100 non-HLA loci were shown to explain 4.7% of the heritability of RA in Asians. The majority of the heritability is still missing.' $input ' have xxxxxx function which might be invovled in  pathology of RA. therefore, In this study, we conducted assocation study to investigate the role of xx and its paralog genes and in RA. in the first stage, we colllected 1078 seropositive RA and 1045 matched control. xxx SNPs in xx, xxx, xxx, xx were genotyped. We found SNPs rsxxx in xxx was signifciantly associated with RA, P=xxx, 95%CI.' >> readme.txt
echo "Run title:" $rs "in" $input "and Seropositive Rheumatoid Arthritis" >> readme.txt
echo "Reference: https://github.com/CNAID/Publication/blob/master/2018/1-s2.0-S155608641830090X-main.pdf" >>readme.txt
echo "Reference: https://github.com/CNAID/Publication/blob/master/2019/nejmoa18015622019.pdf" >>readme.txt
