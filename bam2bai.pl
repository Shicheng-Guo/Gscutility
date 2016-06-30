#!/usr/bin/perl -w

# Transfer Bam to Fastq with samtools command
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-02-19
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
