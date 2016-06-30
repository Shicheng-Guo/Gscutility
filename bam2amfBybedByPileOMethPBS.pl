#!/usr/bin/perl

#bin/perl -w Bam to sort Bam with samtools in TSCC
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-01-19

use strict;
use Cwd;
my $dir=getcwd;
my @file=glob("/home/shg047/oasis/monod/BAM/rename/*.bam");
mkdir "MF_PileOMeth" if ! -e "MF_PileOMeth";
my $ppn=1;
my $walltime="7:00:00";
my $queue="home-k4zhang"; # hotel, home-k4zhang, glean

foreach my $file (@file){
        my @tmp=split/\//,$file;
    my ($sample,undef)=split /.bam/,$tmp[$#tmp];
    my $job_file_name = $sample.".bam2bedGraph.PileOMeth.job";
    open(OUT, ">./MF_PileOMeth/$job_file_name") || die("Error in opening file $job_file_name.\n");
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $job_file_name\n";
    print OUT "#PBS -q $queue\n";
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime\n";
    print OUT "#PBS -o ".$file.".bam2bedGraph.PileOMeth.log\n";
    print OUT "#PBS -e ".$file.".bam2bedGraph.PileOMeth.err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shicheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $dir\n";
    print OUT "PileOMeth extract --keepDupes --keepSingleton --keepDiscordant --mergeContext -p 5 -q 10 --minDepth 5 /home/shg047/oasis/db/hg19/hg19.fa $file -o ./MF_PileOMeth/$sample\n";
    system("qsub ./MF_PileOMeth/$job_file_name");
}
