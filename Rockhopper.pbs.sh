#########################################################################################################
##### RNA-seq data to reveal novel response mechanism to bacterial within host wound tissues ###########
#########################################################################################################
## 02/04/2020
cd ~/hpc/project/RnaseqBacterial/extdata/rnaseq
wget http://cs.wellesley.edu/~btjaden/Rockhopper/download/current/Rockhopper.jar
cd ~/hpc/project/RnaseqBacterial/extdata/rnaseq
genome_DIR1=~/hpc/project/RnaseqBacterial/extdata/rnaseq/Rockhopper_Results/genomes/Staphylococcus_aureus_subsp__aureus_USA300_FPR3757
mkdir temp
mkdir diamond
for i in $(ls *.fastq.gz | rev | cut -c 17- | rev | uniq)
do
echo $i
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=12 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Xmx1200m -cp Rockhopper.jar Rockhopper -g $genome_DIR1 $i\_R1_001.fastq.gz%$i\_R2_001.fastq.gz -o ./Rockhopper_Results/$i >> $i.job
echo pear -f $i\_R1_001.fastq.gz -r $i\_R2_001.fastq.gz -o $i >> $i.job
echo zcat $i\_R1_001.fastq.gz \> $i.fastq >> $i.job
echo zcat $i\_R2_001.fastq.gz \>\> $i.fastq >> $i.job
echo diamond blastx -d ~/hpc/db/nr -q $i.fastq -o ./diamond/$i >> $i.job
qsub  $i.job
done
