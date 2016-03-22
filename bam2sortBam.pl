#!/usr/bin/perl -w

# Bam to sort Bam with samtools in TSCC
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-01-19

use strict;
use Cwd;
die USAGE() if scalar(@ARGV)<2;
my $dir=getcwd;
my @file=glob("$dir/*.bam");
my $destFold=shift @ARGV;
my $submit=shift @ARGV;
mkdir "$destFold" if ! -e "$destFold";

my $ppn=1;
my $walltime="10:00:00";
my $queue="glean"; # hotel

foreach my $file (@file){
    my @tmp=split /\//,$file;
    my ($sample,undef)=split /.bam/,$tmp[$#tmp];  
    my $job_file_name = $file."sort.job";
    my $status_file = $file.".status";
    my $hapInfo_file = $file.".hapInfo.txt";
    my $curr_dir = $dir;
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $job_file_name\n";
    print OUT "#PBS -q $queue\n";  
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime\n";
    print OUT "#PBS -o ".$file.".bam2sortbam.log\n";
    print OUT "#PBS -e ".$file.".bam2sortbam.err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $curr_dir\n";
    print OUT "samtools sort $file -o ../sortbam/$sample.sort.bam\n";
    print OUT "samtools index ../sortbam/$sample.sort.bam\n";	
    if($submit eq 'submit'){
    system("qsub $job_file_name");
    }
}

sub USAGE{
        print "\nUSAGE: perl $0 Bam_Directory Output_Dir submit\n\n";
        print "Example: perl $0 /media/LTS_60T/Dinh/BAM /media/LTS_60T/Shicheng/BAM submit\n";
}
