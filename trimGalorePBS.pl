#!/usr/bin/perl -w
# Trim_galore the single-end and paired-end gz compressed fastq files
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-01-19

use strict;
use Cwd;
my $dir=getcwd;

die USAGE() if scalar(@ARGV)<2;
my $sample_group_file=shift @ARGV;
my $submit=shift @ARGV;

open F,"$sample_group_file";
while(<F>){
    chomp;
    next if /^\s+$/;
    my @read = split /\t/;
    my $id=$read[0];
    my ($sample1,undef)=split /.fastq.gz|_1.fastq.gz|fq.gz|.gz/,$read[0];

    my $ppn=1;
    my $walltime="7:00:00";
    my $queue="glean"; # hotel, glean, pdafm
    
    my $job_file_name = $id . ".job";
    my $status_file = $id.".status";
    mkdir "../fastq_trim" if ! -e "../fastq_trim";
    
    chomp(my $phredcheck=`perl /home/shg047/bin/checkphred.pl $read[0]`);
    my ($phred)=split /\s+/,$phredcheck;
    $phred="--phred$phred";
    
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $job_file_name\n";
    print OUT "#PBS -q $queue\n";  # glean is free
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime\n";
    print OUT "#PBS -o ".$id.".trim.log\n";
    print OUT "#PBS -e ".$id.".trim.err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $dir\n";
    if(scalar(@read) eq 2){
        my ($sample2,undef)=split /.fastq.gz/,$read[1];
        print OUT "trim_galore $phred --fastqc --illumina --paired --clip_R1 3 --clip_R2 3 --three_prime_clip_R1 3 --three_prime_clip_R2 3 $read[0] $read[1] --output_dir ../fastq_trim\n";
    }else{
        print OUT "trim_galore $phred --fastqc --illumina --clip_R1 3 --clip_R2 3 $read[0] --output_dir ../fastq_trim\n";
    }
   close(OUT);
   if($submit eq 'submit'){
   system("qsub $job_file_name");
   }
}

sub USAGE{
        print "\nUSAGE: perl $0 Sample_Group_File submit\n\n";
        print "Please config the Group_of_Sample mannually";
        print '
Group_of_Sample Format:
HOT243.fq.gz
HOT265.fq.gz
HOT273.fq.gz
TBR34prePlsm.fq.gz
TBR36prePlsm.fq.gz
BMT1.read1.fq.gz        BMT1.read2.fq.gz
BMT2.read1.fq.gz        BMT2.read2.fq.gz
BMT3.read1.fq.gz        BMT3.read2.fq.gz
LTP1.read1.fq.gz        LTP1.read2.fq.gz
'
}
