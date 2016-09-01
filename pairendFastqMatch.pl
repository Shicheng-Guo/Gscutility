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
my @fastq=glob("*fastq.gz *fq.gz *fq *fastq *fa");
foreach my $fastq(@fastq){
        my($sample,undef)=split /[_.]/,$fastq;
        push(@{$sam{$sample}},$fastq);
}

foreach my $sam (sort keys %sam){
if(scalar($sam{$sam})>1){
print OUT "$sam{$sam}[0]\t$sam{$sam}[1]\n"	
}else{
print OUT "$sam{$sam}[0]\n";	
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

Format Example:

Indx01_1.fq          Indx01_2.fq
Indx01_1.fastq       Indx01_2.fastq
Indx01_1.fq.gz       Indx01_2.fq.gz
Indx01.read1.fq.gz   Indx01.read2.fq.gz
--------------------------------------------------------------------------------------------------------------------------

Fastq Match Files was saved in current directory as: FastqMatchConfig.txt


';
}
