for i in `ls *.bedgraph`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/epimarker/bedmethyl >> $i.job
echo j=\"${i//.bed.bedgraph/}\" >> $i.job
echo wigToBigWig $i ~/hpc/db/hg38/hg38.chrome.sizes \$j.bw >> $i.job
# qsub $i.job
done
