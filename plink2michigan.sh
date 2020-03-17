#### plink2vcf and Michigan imputation
cd ~/hpc/rheumatology/RA/RA500
mkdir michigan
plink --bfile RA3000.R5 --list-duplicate-vars ids-only suppress-first
plink --bfile RA3000.R5 --alleleACGT --snps-only just-acgt --exclude plink.dupvar --make-bed --out ./michigan/RA3000_R5
cd michigan
mkdir temp
wget https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar -O conform-gt.24May16.cee.jar
for i in {1..23} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=12 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo plink --bfile RA3000_R5 --chr $i --recode vcf-iid --out RA3000_R5.chr$i >> $i.job
echo bcftools view RA3000_R5.chr$i.vcf -Oz -o RA3000_R5.chr$i.vcf.gz >>$i.job
echo tabix -p vcf RA3000_R5.chr$i.vcf.gz >>$i.job
echo # java -jar ./conform-gt.24May16.cee.jar gt=RA2020-B9.chr$i.vcf.gz match=POS chrom=$i ref=~/hpc/db/hg19/beagle/EAS/chr$i.1kg.phase3.v5a.EAS.vcf.gz  out=RA2020-B9.chr$i.beagle >>$i.job
echo # tabix -p vcf RA2020-B9.chr$i.beagle.vcf.gz >>$i.job
sh $i.job &
done

mkdir temp
for i in {1..22} X
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo unzip -P N1fsvXwtY.ZR6Q chr_$i.zip  >> $i.job
qsub $i.job
done


cd ~/hpc/db/hg19/beagle/EUR/
wget https://raw.githubusercontent.com/zhanxw/checkVCF/master/checkVCF.py
samtools faidx hg19.fa
cd /gpfs/home/guosa/hpc/project/pmrp/cytokine/hg18
python checkVCF.py -r hs37d5.fa -o SchrodiTH17_660W.chr9  SchrodiTH17_660W.chr9.vcf.gz
plink --vcf Schrodi_IL23_IL17_combined_RECAL_SNP_INDEL_variants.VA.chr19.Minimac4.vcf.gz --double-id --freq --make-bed --out Schrodi_IL23_IL17.chr19
wget https://faculty.washington.edu/browning/beagle/beagle.16May19.351.jar -O beagle.16May19.351.jar
wget https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar -O conform-gt.24May16.cee.jar
java -jar ./conform-gt.24May16.cee.jar gt=SchrodiTH17_660W.chr22.vcf.gz chrom=22 ref=~/hpc/db/hg19/beagle/EUR/chr22.1kg.phase3.v5a.EUR.vcf.gz  out=SchrodiTH17_660W.chr22.beagle.vcf.gz
#############################################
cd ~/hpc/db/hg19/beagle
for i in {1..22} X Y
do
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr$i.1kg.phase3.v5a.vcf.gz
done
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/20140625_related_individuals.txt
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_male_samples_v3.20130502.ALL.panel
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples.20130502.ALL.ped
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel
mkdir EUR
mkdir EAS

grep EUR integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}'> EUR.List.txt
grep EAS integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > EAS.List.txt

mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
# echo tabix -p vcf chr$i.1kg.phase3.v5a.vcf.gz >> $i.job
# echo bcftools view chr$i.1kg.phase3.v5a.vcf.gz -S EUR.List.txt -Oz -o ./EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz >>$i.job
echo bcftools view chr$i.1kg.phase3.v5a.vcf.gz -S EAS.List.txt -Oz -o ./EAS/chr$i.1kg.phase3.v5a.EAS.vcf.gz >>$i.job
echo bcftools norm -d all -m-both ./EAS/chr$i.1kg.phase3.v5a.dedup.EAS.vcf.gz -Oz -o ./EAS/chr$i.1kg.phase3.v5a.dedup.norm.EAS.vcf.gz  >>$i.job
qsub $i.job
done

mkdir temp
mkdir plink
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo # gunzip chr$i.info.gz >>$i.job
echo # bcftools view -i \'R2\>0.6\|TYPED=1\|TYPED_ONLY=1\' -Oz chr$i.dose.vcf.gz -Oz -o chr$i.dose.filter.vcf.gz >>$i.job
echo # plink --vcf chr$i.dose.vcf.gz --make-bed --out chr$i >>$i.job
echo # plink --bfile chr$i --list-duplicate-vars --out chr$i >>$i.job
echo plink --bfile chr$i --exclude chr$i.dupvar --make-bed --out ./plink/chr$i >> $i.job
qsub $i.job
done

### revise the phentoypes with R script
data1<-read.table("/gpfs/home/guosa/hpc/rheumatology/RA/RA500/RA2020-B9.fam")
file=list.files(pattern="*.fam")
for(i in file){
data2<-read.table(i)
data2[,6]<-data1[match(data2[,1],data1[,2]),6]
out=paste(i,"new",sep=".")
write.table(data2,file=out,sep=" ",quote=F,col.names=F,row.names=F)
cmd= paste("mv",out,i,sep=" ")
print(cmd)
system(cmd)
}

mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo plink --bfile chr$i --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --counts --out chr$i.counts >> $i.job
echo plink --bfile chr$i --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --out chr$i.freq >> $i.job
qsub $i.job
done

wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/high-LD-regions.txt

plink --bfile chr1 --allow-no-sex --maf 0.05 --merge-list merge.txt --make-bed --out RA500

plink --bfile RA500 --exclude high-LD-regions.txt --range --indep-pairwise 50 5 0.2
plink --bfile RA500 --allow-no-sex  --extract plink.prune.in --genome --make-bed --out temp
plink --bfile temp --pca --maf 0.05 --memory 40000 --cluster --mds-plot 4 --out RA500
plink --bfile RA500 --allow-no-sex --logistic --threads 31 --covar RA500.eigenvec --covar-number 1-4 --adjust --out RA500
plink --bfile RA500 --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --counts --out RA500.counts
plink --bfile RA500 --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --out RA500.freq 

for i in {2..22}
do
echo chr$i >> merge.txt
done

cd /gpfs/home/guosa/hpc/rheumatology/RA/RA500/michigan
gunzip chr22.info.gz
plink --vcf chr22.dose.vcf.gz --make-bed --out chr22
plink --bfile chr22 --list-duplicate-vars --out chr22
plink --bfile chr22 --exclude plink.dupvar --make-bed --out ../plink/chr22

plink --bfile RA2020 --mind 0.05 --make-bed --out RA2020-B1
plink --bfile RA2020-B1 --geno 0.95 --make-bed --out RA2020-B2
plink --bfile RA2020-B2 --maf 0.01 --make-bed --out RA2020-B3
plink --bfile RA2020-B3 --hwe 0.00001 --make-bed --out RA2020-B4
plink2 --bfile RA2020-B4 --king-cutoff 0.125
plink2 --bfile RA2020-B4 --remove plink2.king.cutoff.out.id --make-bed -out RA2020-B5
plink --bfile RA2020-B5 --check-sex
plink --bfile RA2020-B5 --impute-sex --make-bed --out RA2020-B6
plink --bfile RA2020-B6 --check-sex
grep PROBLEM plink.sexcheck | awk '{print $1,$2}' > sexcheck.remove
plink --bfile RA2020-B6 --remove sexcheck.remove --make-bed --out RA2020-B7
plink --bfile RA2020-B7 --test-missing midp 
awk '$5<0.000001{print}' plink.missing | awk '{print $2}' > missing.imblance.remove
plink --bfile RA2020-B7 --exclude missing.imblance.remove --make-bed --out RA2020-B8
plink --bfile RA2020-B8 --assoc mperm=5000 --adjust gc --threads 31
plink --bfile RA2020-B8 --pca --threads 31
plink --bfile RA2020-B8 --logistic --covar plink.eigenvec --covar-number 1-20 --adjust

wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200121.zip