cd /gpfs/home/guosa/hpc/methylation/pancrease/medip/fastq
mkdir temp
for i in $(ls *.fastq.gz | rev | cut -c 12- | rev | uniq)
do
echo $i
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=6 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bowtie2 -x /gpfs/home/guosa/hpc/db/hg19/bowtie2/hg19 -1 $i\_1.fastq.gz -2 $i\_2.fastq.gz -S $i.sam >> $i.job
echo samtools view -bS $i.sam \> $i.bam >> $i.job
echo samtools sort $i.bam -o $i.sorted.bam >> $i.job
echo samtools mpileup -uf ~/hpc/db/hg19/hg19.db $i.sorted.bam \| bcftools view -Ov - \> $i.bcf >> $i.job
qsub  $i.job
done
