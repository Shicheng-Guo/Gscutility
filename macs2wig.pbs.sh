cd /gpfs/home/guosa/hpc/methylation/pancrease/medip/bam
cd /home/guosa/hpc/methylation/pancrease/medip/bam
cd /gpfs/home/guosa/hpc/methylation/pancrease/medip/bam/ucsc
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz -O chromInfo.txt.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/bedGraph2wig.pl -O bedGraph2wig.pl
mkdir temp
for i in $(ls *.bam.FE.bdg.wig | rev | cut -c 16- | rev | uniq)
do
echo $i
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=6 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
# echo bowtie2 -p 6 -x /gpfs/home/guosa/hpc/db/hg19/bowtie2/hg19 -1 $i\_1.fastq.gz -2 $i\_2.fastq.gz -S $i.sam >> $i.job
# echo samtools view -bS $i.sam \> $i.bam >> $i.job
# echo samtools sort $i.bam -o $i.sorted.bam >> $i.job
# echo samtools mpileup -uf ~/hpc/db/hg19/hg19.db $i.sorted.bam \| bcftools view -Ov - \> $i.bcf >> $i.job
# echo samtools depth $i.sorted.bam \> $i.wig >> $i.job
# echo macs2 callpeak -t $i -f BAM -g hs -n $i -B -q 0.01 >> $i.job
# echo  macs2 bdgcmp -t $i\_treat_pileup.bdg -c $i\_control_lambda.bdg -o ./ucsc/$i.FE.bdg -m FE >>$i.job
# echo  macs2 bdgcmp -t $i\_treat_pileup.bdg -c $i\_control_lambda.bdg -o ./ucsc/$i.logLR.bdg -m logLR -p 0.00001 >>$i.job
# echo sh bdg2bw.sh ./ucsc/$i.FE.bdg ./ucsc/chromInfo.txt >>$i.job
# echo bigWigToWig ./ucsc/$i.FE.bw ./ucsc/$i.wig >> $i.job
echo perl bedGraph2wig.pl --bedgraph $i.bam.FE.bdg --wig $i.wig --step 150 >> $i.job
echo gzip -c $i.wig \> $i.wig.gz >> $i.job
qsub  $i.job
done
