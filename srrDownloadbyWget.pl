

#!/usr/bin/perl
use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;

my $file=shift @ARGV;
my $submit=shift @ARGV;
open F,$file;
my %list;
while(<F>){
if(/(SRR\d+).*(ftp.*gz)/){
$list{$1}=$2;
}
}

my @job=glob("*.job");
my @fastq=glob("*.gz");

foreach my $job(@job){
my ($sample,undef)=split/\./,$job;
my $exist=1;
foreach my $fastq(@fastq){
$exist=0 if $fastq =~/$sample/;
next;
}
# print "$list{$sample}\n" if $exist eq 1;
if($exist eq 1){
my $id=$sample;
my $curr_dir = $dir;
my $queue="hotel";
my $walltime="168:00:00";
my $ppn=1;
my $job_file_name="$sample.job";
open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
print OUT "#!/bin/csh\n";
print OUT "#PBS -N $id\n";
print OUT "#PBS -q $queue\n";  # glean is free, pdafm
print OUT "#PBS -l nodes=1:ppn=$ppn\n";
print OUT "#PBS -l walltime=$walltime\n";
print OUT "#PBS -o ".$id.".log\n";
print OUT "#PBS -e ".$id.".err\n";
print OUT "#PBS -V\n";
print OUT "#PBS -M shicheng.guo\@gmail.com \n";
print OUT "#PBS -m abe\n";
print OUT "#PBS -A k4zhang-group\n";
print OUT "cd $curr_dir\n";
print OUT "wget $list{$sample}\n";
close(OUT);
if($submit eq 'submit'){
system("qsub $job_file_name");
}
}
}
