for i in `ls *.bam`
do
j=${i/.bam/}
echo \#PBS -N $i.sort  > $i.job
echo \#PBS -l nodes=1:ppn=6 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/bam/pool/bam >> ${i}.job
echo cd $(pwd) >> ${i}.job
echo samtools sort -n -@ 6 $i $j.nsort >> ${i}.job
echo $i.job
qsub $i.job
done
