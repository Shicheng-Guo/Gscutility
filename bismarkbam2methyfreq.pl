#!/usr/bin/perl -w

# extract methBedGraph by Bismark from unsorted bam file (both for single and pair-end. reads should be sort by reads name for pair-end bam)
# Good habit to creat parameter table for a project
# Contact: Shicheng Guo
# Version 1.3
# Update: Feb/28/2016

use strict;
use Cwd;
use Getopt::Long;
my ($input,$submit,$genome,$server,$help,$queue,%ppn,%multicore,%walltime,$is_pairedEnd);
## default parameters
$input ="saminfo.txt";
$submit='nonsubmit';
$genome="hg19";
$server="TSCC";
my $dir=getcwd;
chdir $dir;
my @bamFile=glob("*.bam");
my $curr_dir = $dir;

## parse parameters
GetOptions ( "input=s"   => \$input,          # string
             "genome=s" => \$genome,          # string   
             "server=s" => \$server,          # flag
             "help" => \$help,                # help
	         "queue=s" => \$queue,            # queue(hotel,pdafm,condo)
	    )          
or die("Error in command line arguments\n");
my $BismarkRefereDb;
if($server eq "TSCC"){
	if($genome eq "hg19"){
	$BismarkRefereDb="/home/shg047/db/hg19/bismark/";
	}elsif($genome eq "hg38"){
	$BismarkRefereDb="/home/shg047/db/hg38/bismark/";
	}else{
	warn("print assign genome version (in TSCC)to the script: hg19? hg38? mm9? mm10?");	
	}
}elsif($server eq "GM"){
	if($genome eq "hg19"){
	$BismarkRefereDb="/media/Home_Raid1/shg047/db/hg19/bismark";
	}elsif($genome=="hg38"){
	$BismarkRefereDb="/home/shg047/db/hg38/bismark/";
	}else{
	warn("print assign genome version (in Genome-miner)to the script: hg19? hg38? mm9? mm10?");	
	}
}

    %walltime=(
    "hotel" => "168:00:00",
    "condo" => "8:00:00",
    "pdafm" => "72:00:00",
	"glean" => "168:00:00",
    );
    
    %ppn=(
    hotel => "16",
    pdafm => "32",
    glean => "6",
    condo => "16",
    );

    %multicore=(
    hotel => "6",
    pdafm => "12",
    glean => "1",
    condo => "6",
	glean => "6",
    );
	my $nodes=1;

open F, $input;
while(<F>){
chomp;
my $bamfile=$_;
next if $bamfile=~/.tmp./;

if($bamfile=~/_val/){
$is_pairedEnd="paired";
}else{
$is_pairedEnd="single";
}

my ($sampleid,undef)=split /.bam/,$bamfile;
my $job_file_name = $bamfile.".bam2mf.job";
my $status_file = $bamfile.".status";
mkdir "../methyfreq" if ! -e "../methyfreq";
open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
my($sample,undef)=split /.fq.gz|.fq.gz|.bam|.sort.bam|.sortc.bam|sortn.bam|.fastq.gz/,$bamfile; 	
print OUT "#!/bin/csh\n";
print OUT "#PBS -N $sample\n";
print OUT "#PBS -q $queue\n";  # glean is free
print OUT "#PBS -l nodes=$nodes:ppn=$ppn{$queue}\n";
print OUT "#PBS -l walltime=$walltime{$queue}\n";
print OUT "#PBS -o $sample.log\n";
print OUT "#PBS -e $sample.err\n";
print OUT "#PBS -V\n";
print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
print OUT "#PBS -m abe\n";
print OUT "#PBS -A k4zhang-group\n";
print OUT "cd $curr_dir\n";
if($is_pairedEnd eq 'paired'){                                                      # BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
# print OUT "filter_non_conversion --paired ../bam/$bamfile\n";                       # BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.nonCG_filtered.bam
# print OUT "rm ../bam/$sample\.nonCG_removed_seqs.bam\n";                            # BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.nonCG_removed_seqs.bam
# print OUT "bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 1 --multicore 3 --paired-end --bedGraph --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../bam/$sample.nonCG_filtered.bam\n";
# print OUT "coverage2cytosine --merge_CpG ../methyfreq/ --genome_folder $BismarkRefereDb -o ../methyfreq/$sample\_val_1_bismark_bt2_pe.nonCG_filtered.bismark.CpG_merge.cov";
 print OUT "bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 1 --multicore 6 --paired-end --bedGraph --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../bam/$bamfile\n";
 print OUT "coverage2cytosine --merge_CpG ../methyfreq/ --genome_folder $BismarkRefereDb -o ../methyfreq/$sample\n";
}else{
# print OUT "filter_non_conversion ../bam/$bamfile\n";                                # CTR101_trimmed.fq.gz_bismark_bt2.sort.bam
# print OUT "rm ../bam/$sample\.nonCG_removed_seqs.bam\n";                            # CTR101_trimmed.fq.gz_bismark_bt2.sort.nonCG_removed_seqs.bam
# print OUT "bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 1 --multicore 8 --bedGraph --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../bam/$sample.nonCG_filtered.bam\n";
# print OUT "coverage2cytosine --merge_CpG ../methyfreq/ --genome_folder $BismarkRefereDb -o ../methyfreq/$sample\_val_1_bismark_bt2_pe.nonCG_filtered.bismark.CpG_merge.cov";
 print OUT "bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 1 --multicore 6 --bedGraph --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../bam/$bamfile\n";
 print OUT "coverage2cytosine --merge_CpG ../methyfreq/ --genome_folder $BismarkRefereDb -o ../methyfreq/$sample\n";
}
close OUT;
}

sub USAGE{
print "\nUSAGE: $0 Bam_Directory submit\n";
print "\nUSAGE: $0 ./ submit\n";
print "\nUSAGE: Job Files will be created in BAM_Directory\n";
}
