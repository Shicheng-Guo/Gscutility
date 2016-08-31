#!/usr/bin/perl -w
use strict;
use Cwd;
my $dir=getcwd;
my $submit=shift @ARGV;
my @fastq=glob("*gz");
foreach my $fastq(@fastq){
        my $job_file_name = $fastq.".job";
        open JOB_FILE, ">$job_file_name" || die("Error in opening file $job_file_name.\n");
        print JOB_FILE "#!/bin/csh\n";
        print JOB_FILE "#PBS -q hotel\n";
        print JOB_FILE "#PBS -l nodes=1:ppn=1\n";
        print JOB_FILE "#PBS -l walltime=9:00:00\n";
        print JOB_FILE "#PBS -o ".$fastq.".log\n";
        print JOB_FILE "#PBS -e ".$fastq.".err\n";
        print JOB_FILE "#PBS -V\n";
        print JOB_FILE "#PBS -M shicheng.guo\@gmail.com \n";
        print JOB_FILE "#PBS -m abe\n";
        print JOB_FILE "#PBS -A k4zhang-group\n";
        print JOB_FILE "cd $dir\n";
        print JOB_FILE "gunzip $fastq\n";
        close(JOB_FILE);
        print "$job_file_name created...\n";
        if($submit eq "submit"){
        system("qsub $job_file_name");
        }
}
