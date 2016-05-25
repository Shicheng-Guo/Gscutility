#!/usr/bin/perl -w
use strict;
use Cwd;
# Extract Pileup within specfic genomic regions (CDS, Exon)
# Set Pileup and GenomicInterval regions
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-03-16

my $phred_base=5;
my $scripts_dir="/home/shg047/software/BisReadMapper/src";
my $samtools="/home/shg047/software/samtools-0.1.18/samtools";
my $bamUtils="/home/shg047/software/bamUtil-1.0.12/bin/bam";
my $cpg_list="/home/shg047/oasis/db/hg19_lambda.cpg.positions.txt";
my $ref_fai="/home/shg047/oasis/db/hg19_lambda.fa.fai";
my $ref_fa="/home/shg047/oasis/db/hg19_lambda.fa";
my $submit="NO";
my $bamDir;

die &USAGE if scalar(@ARGV)<2;
$bamDir=shift @ARGV;
$submit=shift @ARGV;

my @bamFile=glob("$bamDir/*bam");
foreach my $sam (@bamFile){
	my @tmp=split/\//,$sam;
	my $id=$tmp[-1];
    my $job_file_name = $id . ".job";
    my $status_file = $id.".status";
    my $hapInfo_file = $id.".hapInfo.txt";
    my $curr_dir = getcwd;
    open(JOB_FILE, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
    print JOB_FILE "#!/bin/csh\n";
    print JOB_FILE "#PBS -q glean\n";
    print JOB_FILE "#PBS -l nodes=1:ppn=1\n";
    print JOB_FILE "#PBS -l walltime=8:00:00\n";
    print JOB_FILE "#PBS -o ".$id.".log\n";
    print JOB_FILE "#PBS -e ".$id.".err\n";
    print JOB_FILE "#PBS -V\n";
    print JOB_FILE "#PBS -M shihcheng.guo\@gmail.com \n";
    print JOB_FILE "#PBS -m abe\n";
    print JOB_FILE "#PBS -A k4zhang-group\n";
    print JOB_FILE "cd $curr_dir\n";
    my $cmd = "$samtools mpileup -BA -f $ref_fa $sam | $scripts_dir/extractMethyl.pl $cpg_list $phred_base > $id.methylFreq";
    print JOB_FILE "$cmd\n";
    close(JOB_FILE);
    print "Job file $job_file_name is created...\n";
    if($submit eq 'submit'){
    system("qsub $job_file_name");
    }
}

sub USAGE{
	print "\nUSAGE: perl $0 bamDir submit\n";
	print 'Enter to destination directory and run the script to extract methylation freqency\n';
	print "# Version 1.3\n";
    print "# Update: 2016-03-16\n";
    print "# Samtools Version: samtools-0.1.18\n";
    print "# $bamUtils Version: bamUtil-1.0.12\n";
    
}
