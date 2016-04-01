#!/usr/bin/perl -w

# Bam to sort Bam with samtools in TSCC
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-01-19

use strict;
use Cwd;
my $dir=getcwd;
my $submit =shift @ARGV;
my @file=glob("SRS*.txt");
my $ppn=1;
my $walltime="14:00:00";
my $queue="glean"; # hotel

foreach my $file (@file){
    my ($sample,undef)=split /.txt/,$file;
    my $job_file_name = $file.".job";
    my $curr_dir = $dir;
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $job_file_name\n";
    print OUT "#PBS -q $queue\n";
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime\n";
    print OUT "#PBS -o ".$file.".SRR2SRS.log\n";
    print OUT "#PBS -e ".$file.".SRR2SRS.err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shicheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $curr_dir\n";
    print OUT "samtools merge -h inh.sam -n -b $file $sample.nsort.bam\n";
    print OUT "samtools sort -o $sample.csort.bam $sample.nsort.bam\n";
#    print OUT "samtools index $sample.csort.bam\n";
#    print OUT "bismark_methylation_extractor --paired-end --bedGraph --multicore 3 --ignore 3 --ignore_3prime 3 --ignore_r2 5 --ignore_3prime_r2 5 --gzip --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  $sample.nsort.bam\n";
#    print OUT ""
#    print OUT "perl ~/bin/bam2hapInfo2PBS.pl ../  $sample.nsort.bam\n";
     if ($submit eq 'submit'){
     system("qsub $job_file_name");
        }
}
