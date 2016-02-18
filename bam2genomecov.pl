 #!/bin/csh
 #PBS -N GenomeCovDepth
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=16
 #PBS -l walltime=72:00:00
 #PBS -o CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.log
 #PBS -e CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.err
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/Haib/sortBam
 # samtools cat -h header.sam -o haib.merge.bam *sort.bam
 samtools sort -@ 16 haib.encode.merge.bam -o haib.merge.sort.bam
 samtools index haib.merge.sort.bam
 bedtools genomecov -bg -split -ibam haib.merge.sort.bam >   haib.merge.bam.pool.bed
 awk '$4>9 { print $1"\t"$2"\t"$3}'  haib.merge.bam.pool.bed | bedtools merge -d 10 -i - > haib.RD10.genomecov.bed
 awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' haib.RD10.genomecov.bed > haib.RD10_80up.genomecov.bed
