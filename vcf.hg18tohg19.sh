### plink_hg18tohg19
### vcf_hg18tohg19
plink="Th17Set1"
wget https://faculty.washington.edu/browning/beagle/beagle.16May19.351.jar -O beagle.16May19.351.jar
wget https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar -O conform-gt.24May16.cee.jar
mkdir temp
for i in {1..23} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=12 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo plink --bfile $plink --chr $i --recode vcf-iid --out $plink.chr$i >> $i.job
echo java -jar ./conform-gt.24May16.cee.jar gt=$plink.chr$i.vcf chrom=$i ref=~/hpc/db/hg19/beagle/EUR/chr$i.1kg.phase3.v5a.dedup.norm.EUR.vcf.gz out=$plink.chr$i.beagle >> $i.job
# echo bcftools norm -m-any chr$i.1kg.phase3.v5a.EUR.vcf.gz -Oz -o chr$i.1kg.phase3.v5a.dedup.EUR.vcf.gz  >>$i.job
# echo bcftools norm -d all chr$i.1kg.phase3.v5a.dedup.EUR.vcf.gz -Oz -o chr$i.1kg.phase3.v5a.dedup.norm.EUR.vcf.gz  >>$i.job
# echo plink --bfile S_Hebbring_Unr.Guo --recode vcf --chr $i --snps-only just-acgt --out ./beagle/S_Hebbring_Unr.Guo.Forward.chr$i >> chr$i.job
# echo unzip -P fptULXgh3wF9Q chr_$i.zip  >> $i.job
qsub $i.job
done
