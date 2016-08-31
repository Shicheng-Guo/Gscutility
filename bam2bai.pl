#!/usr/bin/perl -w

# Transfer Bam to Fastq with samtools command
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-08-19
use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;
my @file=glob("*.bam");
foreach my $file(@file){
 my ($sample,undef)=split /\./,$file;
 open OUT,">$file.bam2bai.job";
 print OUT "#!/bin/csh\n";
 print OUT "#PBS -n $file.bam2bai.job\n";
 print OUT "#PBS -q glean\n";  # glean,condo,hotel
 print OUT "#PBS -l nodes=1:ppn=1\n";
 print OUT "#PBS -l walltime=7:00:00\n";
 print OUT "#PBS -o ".$file.".bam2bai.log\n";
 print OUT "#PBS -e ".$file.".bam2bai.err\n";
 print OUT "#PBS -V\n";
 print OUT "#PBS -M shicheng.guo\@gmail.com \n";
 print OUT "#PBS -m abe\n";
 print OUT "#PBS -A k4zhang-group\n";
 print OUT "cd $dir\n";
 print OUT "samtools index $file";
 close OUT;
 system("qsub $file.bam2bai.job");
}

sub USAGE{
print "\nUSAGE: perl $0  submit\n";
print '
--------------------------------------------------------------------------------------------------------------
This command will index for sorted bam files under PBS system. 
Example:  perl bam2bai.pl submit

Previous Script: 
1, Download or Prepare SRA Download Configure  (http://www.ebi.ac.uk/ena/data/view/SRP028600)
2, perl fastqDownload.pl SamConfig.txt 8 submit   (fastq-dump --split-files --gzip)
3, trim_galore --phred33 --fastqc --stringency 3 -q 20 --trim1 --length 20 --gzip --clip_R1 2 --three_prime_clip_R1 2 --illumina SRR299055_1.fastq.gz --output_dir ../fastq_trim
4, bismark_genome_preparation ./     # merge all the fa to one file (mm9.fa or hg19.fa)
5, bismark --bowtie2 --phred33-quals --fastq -L 30 -N 1 --multicore 2 /home/shg047/oasis/db/mm9 -1 ../fastq_trim/SRR610033_1_val_1.fq.gz -2 ../fastq_trim/SRR610033_2_val_2.fq.gz  -o ../bam	
--------------------------------------------------------------------------------------------------------------
';


}
