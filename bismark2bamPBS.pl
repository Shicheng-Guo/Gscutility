#!/usr/bin/perl -w

# bismark to alignment single and pair-end fastq in same project
# Good habit to creat parameter table for a project
# Contact: Shicheng Guo
# Version 1.3
# Update: Jan/28/2016

use strict;
use Cwd;
my $dir=getcwd;
die USAGE() if scalar(@ARGV<2);

my $sample_group_file=shift @ARGV;
my $submit=shift @ARGV;

open F,"$sample_group_file";
while(<F>){
    chomp;
    next if /^\s+$/;
    my @read = split /\t/;
    my $id=$read[0];
    my ($sample1,undef)=split /.fastq.gz/,$read[0];
    my $project="DE";
    my $analysis="AL";
    my $job_file_name = $id . ".job";
    my $status_file = $id.".status";
    my $curr_dir = $dir;
    
    my $ppn=16;
    my $multicore=6;
    my $walltime="167:00:00";
    my $queue="hotel"; # hotel,pdafm,condo
    # chomp(my $phredcheck=`perl /home/shg047/bin/checkphred.pl $read[0]`);
    # my ($phred)=split /\s+/,$phredcheck;
    my $phred=33;
    $phred="--phred$phred";
    
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");   
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $job_file_name\n";
    print OUT "#PBS -q $queue\n";  # glean is free
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime\n";
    print OUT "#PBS -o ".$id.".log\n";
    print OUT "#PBS -e ".$id.".err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shicheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $curr_dir\n";
    print "$job_file_name\n";
    if(scalar(@read) eq 2){
    my ($sample2,undef)=split /.fq.gz/,$read[1];
    # print OUT "bismark --bowtie2 $phred-quals --fastq -L 30 -N 1 --multicore 6 /home/shg047/db/hg19/meth/bismark -1 ../fastq_trim/$sample1\_val_1.fq.gz -2 ../fastq_trim/$sample2\_val_2.fq.gz  -o ../bam";
    print OUT "bismark --bowtie2 $phred-quals --fastq -L 30 -N 1 --multicore 6 /home/shg047/db/mm9 -1 ../fastq_trim/$sample1\_val_1.fq.gz -2 ../fastq_trim/$sample2\_val_2.fq.gz  -o ../bam";
    }else{
    # print OUT "bismark --bowtie2 $phred-quals --fastq -L 30 -N 1 --multicore 6 /home/shg047/db/hg19/meth/bismark ../fastq_trim/$sample1\_trimmed.fq.gz -o ../bam2";  
    print OUT "bismark --bowtie2 $phred-quals --fastq -L 30 -N 1 --multicore 6 /home/shg047/db/mm9 ../fastq_trim/$sample1\_trimmed.fq.gz -o ../bam";  
    }
   close(OUT);
   if($submit eq 'submit'){
   system("qsub $job_file_name");
   }
}

sub USAGE{
print "\nUSAGE: $0 Sam_Config_File submit\n";
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
';
}
