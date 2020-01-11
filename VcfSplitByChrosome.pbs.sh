mkdir temp
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/contigReplace.pl -O contigReplace.pl
for i in {1..23} X
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view dbSNP153.GRCh38p12b.vcf.gz -r $i -Ov -o ./chr/dbSNP153.chr$i.hg38.vcf >>$i.job
echo tabix -p vcf chr$i.dose.contig.vcf.gz >>$i.job
qsub $i.job
done
