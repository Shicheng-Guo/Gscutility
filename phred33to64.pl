#!/usr/bin/perl -w
use strict;
use Cwd;
my $input=shift @ARGV;
my $dir=getcwd;

my @fastq=glob("*fastq|*fq");
foreach my $fastq(@fastq){
	my $id=split /\./,$fastq;
	my $job_file_name = $id . ".job";
    my $status_file = $id.".status";
    
    chomp(my $phredcheck=`perl /home/shg047/bin/checkphred.pl $fastq`);
    my ($phred)=split /\s+/,$phredcheck;    
    if($phred eq 33){
    	open JOB_FILE, ">$job_file_name" || die("Error in opening file $job_file_name.\n");
        print JOB_FILE "#!/bin/csh\n";
        print JOB_FILE "#PBS -q hotel\n";
        print JOB_FILE "#PBS -l nodes=1:ppn=1\n";
        print JOB_FILE "#PBS -l walltime=12:00:00\n";
        print JOB_FILE "#PBS -o ".$id.".log\n";
        print JOB_FILE "#PBS -e ".$id.".err\n";
        print JOB_FILE "#PBS -V\n";
        print JOB_FILE "#PBS -M shihcheng.guo\@gmail.com \n";
        print JOB_FILE "#PBS -m abe\n";
        print JOB_FILE "#PBS -A k4zhang-group\n";
        print JOB_FILE "cd $dir\n";
        print JOB_FILE "reformat.sh in=$id out=$id.phred64.fq qin=33 qout=64\n";
        close(JOB_FILE);
        system("qsub $job_file_name");
    }
}
