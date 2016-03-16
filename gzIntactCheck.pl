#! /usr/bin/perl
use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;
my @file=glob("*gz");

foreach my $file(@file){
 my $pbs_job_submit="$file.gzcheck.job";
 open OUT,">$pbs_job_submit";
 print OUT "#!/bin/csh\n";
 print OUT "#PBS -N $file.job\n";
 print OUT "#PBS -q glean\n";  # glean,condo,hotel
 print OUT "#PBS -l nodes=1:ppn=1\n";
 print OUT "#PBS -l walltime=1:00:00\n";
 print OUT "#PBS -o ".$file.".gzcheck.log\n";
 print OUT "#PBS -e ".$file.".gzcheck.err\n";
 print OUT "#PBS -V\n";
 print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
 print OUT "#PBS -m abe\n";
 print OUT "#PBS -A k4zhang-group\n";
 print OUT "cd $dir\n";
 print OUT "gunzip -c $file | head 2> $file.intact\n";
 close OUT;
 system("qsub $pbs_job_submit");
}
