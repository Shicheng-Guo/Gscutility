#!/usr/bin/perl -w

# extract methBedGraph from unsorted bam file (both for single and pair-end. reads should be sort by reads name for pair-end bam)
# Good habit to creat parameter table for a project
# Contact: Shicheng Guo
# Version 1.3
# Update: Feb/28/2016

use strict;
use Cwd;
die "Usage: perl $0 BAM_Directory submit\n" if scalar(@ARGV)<2;
my $dir=shift @ARGV;
my $submit=shift @ARGV;
chdir $dir;
my @bamFile=glob("*.bam");
my $curr_dir = $dir;
my $queue="glean";
my $ppn=6;
my $walltime="7:00:00";

foreach my $bamfile(@bamFile){
my ($sampleid,undef)=split /.bam/,$bamfile;
my $job_file_name = $bamfile."bam2mf.job";
my $status_file = $bamfile.".status";
mkdir "../methyfreq" if ! -e "../methyfreq";
print "$sampleid.sort.bam\n";
open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
print OUT "#!/bin/csh\n";
print OUT "#PBS -N $bamfile\n";
print OUT "#PBS -q $queue\n";
print OUT "#PBS -l nodes=1:ppn=$ppn\n";
print OUT "#PBS -l walltime=$walltime\n";
print OUT "#PBS -o ".$bamfile.".log\n";
print OUT "#PBS -e ".$bamfile.".bam2mf.err\n";
print OUT "#PBS -V\n";
print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
print OUT "#PBS -m abe\n";
print OUT "#PBS -A k4zhang-group\n";
print OUT "cd $curr_dir\n";
my $is_pairedEnd= & is_bamPE("$bamfile");
print "$bamfile: $is_pairedEnd ...\n";
print OUT "bismark_methylation_extractor $is_pairedEnd --bedGraph --multicore 3 --ignore 3 --ignore_3prime 3 --ignore_r2 5 --ignore_3prime_r2 5 --gzip --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  $bamfile\n";
close OUT;
if($submit eq "submit"){
system("qsub $job_file_name");
}
}

sub is_bamPE($){
  my($bamfile)=@_;
  my $is_pairedEnd;
  my $line=`samtools view $bamfile | head -n 1| awk '{print \$2}'`;
  if($line % 2){
  $is_pairedEnd="--paired-end"
  }else{
  $is_pairedEnd="--single-end"
  }
  return $is_pairedEnd;
}


sub USAGE{
print "\nUSAGE: Note: It only works for BAMs Aligned By Bismark!\n"
print "\nUSAGE: $0 Bam_Directory submit\n";
print "\nUSAGE: $0 ./ submit\n"
print "\nUSAGE: Job Files will be created in BAM_Directory\n"
}
