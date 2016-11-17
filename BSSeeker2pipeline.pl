#!/usr/bin/perl -w

# Convert bed3 to bed4
# Contact: Shicheng Guo
# Version 1.3
# Update: Jan/19/2016

use strict;
use Cwd;
my $dir=getcwd;

my $bsseekdir="~/software/BSseeker2";

# build directory
mkdir "../bedgraph" if ! -e "../bedgraph";
mkdir "../bam" if ! -e "../bam";

# DMR => extend 200bp => bedtools getfasta => build bsseeker index
my $targetfa="/media/Home_Raid1/shg047/NAS3/Minghua2016/fa/target.fa";
my $indexdb="/media/Home_Raid1/shg047/software/BSseeker2/bs_utils/reference_genomes/target.fa_bowtie2";


my @file=glob("*fastq");
my $filenum=$#file+1;
print "$filenum files were found in the folder\n";

## pair-end file regroup
my %sample;
foreach my $file(@file){
        my($sample,undef)=split/_/,$file;  # change:separtor to get sample id
    $sample{$sample}=$sample;
}

## Aligment and methylation calling
foreach my $sample(sort keys %sample){
        my $read1="$sample\_R1.fastq";      # change: sample read1
        my $read2="$sample\_R2.fastq";      # change: sample read2
        my $out1=system("python $bsseekdir/bs_seeker2-align.py -1 $read1 -2 $read2 --aligner=bowtie2 -o ../bam/$sample.bam -f bam -g $targetfa");
        my $out2=system("python $bsseekdir/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i ../bam/$sample.bam -o ../bedgraph/$sample --txt --db $indexdb");
        print "python $bsseekdir/bs_seeker2-align.py -1 $read1 -2 $read2 --aligner=bowtie2 -o ../bam/$sample.bam -f bam -g $targetfa\n";
        print "python $bsseekdir/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i ../bam/$sample.bam -o ../bedgraph/$sample --txt --db $indexdb\n\n";
}


