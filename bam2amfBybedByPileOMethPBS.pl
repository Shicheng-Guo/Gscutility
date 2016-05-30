#!/usr/bin/perl -w
use strict;
use Cwd;
my $cur_dir=getcwd;
chdir $cur_dir;

my @bed=glob("/home/shg047/oasis/monod/mhb/split/x*");
my $submit=shift @ARGV;
foreach my $bed (@bed){
    my @tmp=split/\//,$bed;
    my $id=$tmp[-1];
    next if $id=~/LambdaNEB/i;
    my $job_file_name = $id . ".job";
    my $status_file = $id.".status";
        my $hapInfo_file = $id.".hapInfo.txt";
        my $curr_dir = `pwd`;
        open(JOB_FILE, "> $job_file_name") || die("Error in opening file $job_file_name.\n");
        print JOB_FILE "#!/bin/csh\n";
        print JOB_FILE "#PBS -N $id\n";
        print JOB_FILE "#PBS -q glean\n";
        print JOB_FILE "#PBS -l nodes=1:ppn=1\n";
        print JOB_FILE "#PBS -l walltime=14:00:00\n";
        print JOB_FILE "#PBS -o ".$id.".log\n";
        print JOB_FILE "#PBS -e ".$id.".err\n";
        print JOB_FILE "#PBS -V\n";
        print JOB_FILE "#PBS -M shicheng.guo\@gmail.com \n";
        print JOB_FILE "#PBS -m abe\n";
        print JOB_FILE "#PBS -A k4zhang-group\n";
        print JOB_FILE "cd $curr_dir\n";
        my $cmd = "perl PileOMethByBed.pl $bed $id";
        print JOB_FILE "$cmd\n";
        close(JOB_FILE);
        print "Job file $job_file_name created.\n";
        if($submit eq 'submit'){
        system("qsub $job_file_name");
        }
}
