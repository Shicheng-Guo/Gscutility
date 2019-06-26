mkdir ~/hpc/tools/bcftools-1.9/db
cp ~/hpc/db/hg19/hg19.fa ~/hpc/tools/bcftools-1.9/db
perl -p -i -e 's/>chr/>/' ~/hpc/tools/bcftools-1.9/db/hg19.fa
perl -p -i -e 's/>X/>23/' ~/hpc/tools/bcftools-1.9/db/hg19.fa
perl -p -i -e 's/>Y/>24/' ~/hpc/tools/bcftools-1.9/db/hg19.fa
samtools faidx ~/hpc/tools/bcftools-1.9/db/hg19.fa
# REF/ALT total/modified/added:   2601741/397967/114
bcftools norm -t "^24,25,26" -m-any --check-ref s -f ~/hpc/tools/bcftools-1.9/db/hg19.fa All_samples_Exome_QC.clean.norm.vcf.gz -Ov | bcftools annotate -x ID,INFO,FORMAT -I +'%CHROM:%POS' -Oz -o All_samples_Exome_QC.clean.vcf.gz

mkdir chr
mkdir temp
wget https://faculty.washington.edu/browning/beagle/beagle.16May19.351.jar -O beagle.16May19.351.jar
wget https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar -O conform-gt.24May16.cee.jar

inputvcf="All_samples_Exome_QC.clean.vcf.gz"
tabix -f -p vcf All_samples_Exome_QC.clean.vcf.gz
for i in {1..23}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=24 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo mkdir ./temp/chr$i >> $i.job
echo bcftools view -e \'ALT =\"-\" \| REF =\"-\"\' -t $i $inputvcf -Oz -o ./chr/All_samples_Exome_QC.clean.chr$i.vcf.gz >>$i.job
echo java -Djava.io.tmpdir=./temp/chr$i -Xmx64g -jar beagle.16May19.351.jar gt=./chr/All_samples_Exome_QC.clean.chr$i.vcf.gz ref=~/hpc/db/hg19/beagle/EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz map=~/hpc/db/hg19/beagle/plink.chr$i.GRCh37.map out=./chr/All_samples_Exome_QC.chr$i.vcf >>$i.job
qsub $i.job
done
