mkdir ../methyfreq
option1=$(echo --no_overlap --merge_non_CpG --cutoff 1 --multicore 5 --paired-end)
option2=$(echo --bedGraph --ignore 1 --buffer_size 4G --comprehensive)
for i in `ls *bam`
do
j=${i/_L001_R1_001_00_bismark_bt2_pe.bam/}
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/nash/bam >> ${i}.job
echo bismark_methylation_extractor ${option1} ${option2} --output ../methyfreq  ./$i >> ${i}.job
echo $j
qsub $i.job
done
