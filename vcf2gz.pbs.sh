cd /gpfs/home/guosa/hpc/db/hg19/1000Genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >> chr$i.job
echo bgzip -c chr$i.uni.vcf  \> chr$i.vcf.gz >> chr$i.job
echo tabix -p vcf chr$i.vcf.gz >> chr$i.job
qsub chr$i.job
done
