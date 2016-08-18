#!/usr/bin/perl -w

# Collect Single-end and Pired-end Fastq files (Downloaded from SRA database)
# Then Bismark could recongize them to decide single-end or pair-end aligment. 
# Run the script to the Fastq directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-08-19

use strict;
use Cwd;

USAGE();

chdir getcwd;
open OUT,">FastqMatchConfig.txt";
my %sam;
my @fastq=glob("*fastq.gz");
foreach my $fastq(@fastq){
        my($sample,undef)=split /[_.]/,$fastq;
        $sam{$sample}++;
}

foreach my $sam (sort keys %sam){
if($sam{$sam}>1){
print OUT "$sam\_1.fastq.gz\t$sam\_2.fastq.gz\n"	
}else{
print OUT "$sam\_1.fastq.gz\n";	
}
}
close OUT;

sub USAGE{
print "\nUSAGE: perl $0\n";
print '

--------------------------------------------------------------------------------------------------------------------------
This command creat input(configure) files for bismarker alignment. 
Single-end Fastq contains one file each line. 
Pair-end Fastq file contains two files each line.
bismark2bamPBS.pl could recognize this file and take it (FastqMatchConfig.txt) as input to make alignment for BS-seq data.
--------------------------------------------------------------------------------------------------------------------------

Fastq Match Files was saved in current directory as: FastqMatchConfig.txt

';
}
