#!/usr/bin/perl -w
# SRA to FASTQ
# Run the script to the SRA directory
# Contact: Shicheng Guo
# Version 1.3
# Update: OCT/19/2016

use strict;
use Cwd;
my $dir=getcwd;
die "Usage: perl $0 <submit|nosubmit>"; 
my $project="SRA2Fastq";
my $analysis="";
my $ppn=1;
my $walltime="168:00:00";
my $queue="hotel"; 
my $submit=shift @ARGV;

my @SRA=glob("*sra");
foreach my $SRA(@SRA){
    next if /SRRID/;
    next if /^\s+$/;
    my $id=$SRA;
    
    print "$id\n";
    my $job_file_name = $id . ".fastq.download.job";
    my $status_file = $id.".status";
    my $curr_dir = $dir;
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $id\n";
    print OUT "#PBS -q $queue\n";  # glean is free, pdafm
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime\n";
    print OUT "#PBS -o ".$id.".download.log\n";
    print OUT "#PBS -e ".$id.".download.err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shicheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $curr_dir\n";
    print OUT "fastq-dump --split-files --gzip $id\n";    
    close(OUT);
    if($submit eq 'submit'){
    system("qsub $job_file_name");
   }
}
