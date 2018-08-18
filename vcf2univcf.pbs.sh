for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo bcftools norm -D --threads=2 ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  -o chr$i.vcf.gz >> chr$i.job
qsub chr$i.job
done
