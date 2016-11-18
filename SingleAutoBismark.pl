#!/usr/bin/perl -w

# BismarkAlignment.pl: a perl script to run bismark pipeline one by one.
# Run the script trim_galore, adaptor remove, bismark alignment and methylfreq and so on
# Contact: Shicheng Guo
# Version 1.3

use strict;
use warnings;
use Cwd;

my $submit=shift @ARGV;
my $dir=getcwd;
chdir $dir;
my @file=glob("*.fastq.gz");
foreach my $file(@file){
my ($sample,undef)=split /\.fastq/,$file;
open OUT,">$file.bismark.pbs";
chomp(my $phredcheck=`perl /home/shg047/bin/checkphred.pl $file`);

my ($phred)=split /\s+/,$phredcheck;

if(! defined $phred){
print "phred cannot be obtained and will be set to 33\n";
$phred=33;
}
mkdir "../fastq_trim" if ! -e "../fastq_trim";
mkdir "../bam" if ! -e "../bam";
mkdir "../bedgraph" if ! -e "../bedgraph";
mkdir "../sortbam" if ! -e "../sortbam";

$phred="--phred$phred";
print OUT "#!/bin/csh\n";
print OUT "#PBS -q condo\n";
print OUT "#PBS -l nodes=1:ppn=8\n";
print OUT "#PBS -l walltime=8:00:00\n";
print OUT "#PBS -o ".$sample.".log\n";
print OUT "#PBS -e ".$sample.".err\n";
print OUT "#PBS -V\n";
print OUT "#PBS -M shicheng.guo\@gmail.com \n";
print OUT "#PBS -m abe\n";
print OUT "#PBS -A k4zhang-group\n";
print OUT "cd $dir\n";
print OUT "START=\$(date +%s)\n";
print OUT "trim_galore $phred --fastqc --illumina --rrbs $file --output_dir ../fastq_trim\n";
print OUT "bismark --bowtie2 $phred-quals --multicore=2 --fastq -L 32 -N 0 -D 5 -R 1 /home/shg047/oasis/db/hg19/align/bismark ../fastq_trim/$sample\_trimmed.fq.gz -o ../bam\n";
print OUT "bismark_methylation_extractor --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../bedgraph  ../bam/$sample\_trimmed.fq_bismark_bt2.bam";
print OUT "samtools sort ../bam/$sample\_trimmed.fq_bismark_bt2.bam -o ../sortbam/$sample\_trimmed.fq_bismark_bt2.sort.bam\n";
print OUT "samtools index ../sortbam/$sample\_trimmed.fq_bismark_bt2.sort.bam\n";
print OUT "END=\$(date +%s)\n";
print OUT "DIFF=\$(echo \"\$END - \$START\" \| bc)\n";
print OUT "echo \"It takes DIFF=\$DIFF seconds to complete this task...\"\n";


if ($submit eq "submit"){
system("qsub $file.bismark.pbs");
}
}

