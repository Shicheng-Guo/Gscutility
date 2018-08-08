
for i in `ls *.bed`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/epimarker/bedmethyl/test >> $i.job
echo perl bedMethyl2bedgraph.pl $i >> $i.job
echo wigToBigWig $i.bedgraph ~/hpc/db/hg38/hg38.chrom.sizes $i.bw >> $i.job
qsub $i.job
done
