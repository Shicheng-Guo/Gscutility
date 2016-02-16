#!/usr/bin/perl -w

# bismark to alignment single-end fastq in same project
# Good habit to creat parameter table for a project
# Contact: Shicheng Guo
# Version 1.3
# Update: Feb/5/2016

use strict;
use Cwd;
my $dir=getcwd;
die "non-submit: perl bam2sortbam.pl ../Sample.info.txt 0 \nsubmit: perl bam2sortbam.pl ../Sample.info.txt 1\n" if scalar(@ARGV)<1;
my $submit=shift @ARGV;

my @bamFile=glob("*sort.bam");
my $project="Ho";
my $analysis="Lung";
my $curr_dir = $dir;
my $queue="glean";
my $ppn=1;
my $walltime="72:00:00";

foreach my $bamfile(@bamFile){
my ($sampleid,undef)=split /\.sort.bam/,$bamfile;
# next if -e "$sampleid.sort.bam";
my $job_file_name = $bamfile.".methylFreq.job";
my $status_file = $bamfile.".status";

print "$sampleid.sort.bam\n";
open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
print OUT "#!/bin/csh\n";
print OUT "#PBS -n $bamfile\n";
print OUT "#PBS -q $queue\n";
print OUT "#PBS -l nodes=1:ppn=$ppn\n";
print OUT "#PBS -l walltime=$walltime\n";
print OUT "#PBS -o ".$bamfile.".log\n";
print OUT "#PBS -e ".$bamfile.".err\n";
print OUT "#PBS -V\n";
print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
print OUT "#PBS -m abe\n";
print OUT "#PBS -A k4zhang-group\n";
print OUT "cd $curr_dir\n";
# print OUT "samtools sort $bamfile -o $sampleid.sort.bam\n";
# print OUT "samtools index $sampleid.sort.bam\n";
my $is_pairedEnd= & is_bamPE("$sampleid.sort.bam");
print OUT "bismark_methylation_extractor $is_pairedEnd --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  $sampleid.sort.bam";
close OUT;
if($submit){
    system("qsub $job_file_name");
    }
}


sub is_bamPE($){
  my($bamfile)=@_;
  my $is_pairedEnd;
  chomp( my $line=`samtools view $bamfile | head -n 1| awk '{print \$2}'`);
  my $remainder=$line%2;
  print "$line\t$remainder\n";
  if($remainder){
  $is_pairedEnd="--paired-end"
  }else{
  $is_pairedEnd="--single-end"
  }
  return $is_pairedEnd;
}
