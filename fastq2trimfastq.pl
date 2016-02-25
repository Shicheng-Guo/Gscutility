#!/usr/bin/perl -w

# bismark to alignment single-end fastq in same project
# Good habit to creat parameter table for a project
# Contact: Shicheng Guo
# Version 1.3
# Update: Jan/28/2016

use strict;
use Cwd;
my $dir=getcwd;
my $samconfig=shift @ARGV;
my $submit=shift @ARGV;


open F,$samconfig;
while(<F>){
    chomp;
    next if /^\s+$/;
    next if /SRRID/;
    my @read = split /\t/;
    my $fastq1=$read[0];
    my $fastq2=$read[1];
    my ($fastq,undef)=split /\_/,$fastq2;
    my $trimFastq=$read[2];

    my $project="Ho";
    my $analysis="Lung";
    my $job_file_name = $fastq . ".trim.job";
    my $status_file = $fastq.".trim.status";
    my $curr_dir = $dir;

    my $ppn=1;
    my $multicore=1;
    my $walltime="7:00:00";
    my $queue="glean"; # hotel,pdafm,condo
    chomp(my $phredcheck=`perl /home/shg047/bin/checkphred.pl $fastq1`);
    my ($phred)=split /\s+/,$phredcheck;
    $phred="--phred$phred";

    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $fastq\n";
    print OUT "#PBS -q $queue\n";  # glean is free
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime\n";
    print OUT "#PBS -o ".$fastq.".log\n";
    print OUT "#PBS -e ".$fastq.".err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $curr_dir\n";
    mkdir "../fastq_trim" if(! -e "../fastq_trim");
    if(defined $phred){
    	if(! -e "../fastq_trim/$trimFastq"){
        	print OUT "trim_galore $phred --fastqc --stringency 3 -q 20 --trim1 --length 20 --gzip --paired --clip_R1 2 --clip_R2 2 --three_prime_clip_R1 2 --three_prime_clip_R2 2 --illumina $fastq1 $fastq2 --output_dir ../fastq_trim\n";
        	close(OUT);
        	if($submit){
        	system("qsub $job_file_name");
        	}
    	}
 	}
}
