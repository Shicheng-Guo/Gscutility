#!/usr/bin/perl -w

# bismark to alignment single and pair-end fastq in same project
# Good habit to creat parameter table for a project
# Contact: Shicheng Guo
# Version 1.3
# Update: Nov/23/2016

use strict;
use Cwd;
use Getopt::Long;

my $input ="saminfo.txt";
my $submit='nonsubmit';
my $genome="hg19";
my $server="TSCC";

GetOptions ( "input=s"   => \$input,           # string
             "submit=s"  => \$submit,          # flag
             "genome=s" => \$genome,          # string   
             "server=s" => \$server)          # flag
or die("Error in command line arguments\n");

my $dir=getcwd;

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

mkdir "../fastq_trim" if ! -e "../fastq_trim";
mkdir "../bam" if ! -e "../bam";
mkdir "../bedgraph" if ! -e "../bedgraph";
mkdir "../sortbam" if ! -e "../sortbam";
mkdir "../methyfreq" if ! -e "../methyfreq";


open F,$input;
while(<F>){
    chomp;
    next if /^\s+$/;
    my @read = split /\t/;
    
    my $phred;
    if($server eq "GM"){
    chomp(my $phredcheck=`perl /media/Home_Raid1/shg047/bin/checkphred.pl $read[0]`);
    ($phred)=split /\s+/,$phredcheck; 
    $phred=33 if ! defined $phred;
    print "$read[0]\tphred= $phred\n";
    }elsif($server eq "TSCC"){
    chomp(my $phredcheck=`perl ~/bin/checkphred.pl $read[0]`);
    ($phred)=split /\s+/,$phredcheck; 
    $phred=33 if ! defined $phred; 
    print "$read[0]\tphred= $phred\n";
    }
    
    my $nodes=1;
    my $ppn=8;
    my $multicore=2;
    my $curr_dir = $dir;
    my $queue="hotel"; # hotel,pdafm,condo
    my %walltime=(
    "hotel" => "168:00:00",
    "condo" => "8:00:00",
    "pdafm" => "72:00:00",
    );
    
	
    if(scalar(@read) eq 2){
   
    my($sample,undef)=split /_1.fastq.gz|_R1.fastq.gz/,$read[0]; 	
    my($sample1,undef)=split /.fastq.gz/,$read[0];
    my($sample2,undef)=split /.fastq.gz/,$read[1];
    my $job_file_name = "$sample.pbs";
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");   
    
	print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $sample\n";
    print OUT "#PBS -q $queue\n";  # glean is free
    print OUT "#PBS -l nodes=$nodes:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime{$queue}\n";
    print OUT "#PBS -o $sample.log\n";
    print OUT "#PBS -e $sample.err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";

    print OUT "cd $curr_dir\n";    
    print OUT "# fastq-dump --split-files --gzip $sample\n";
	print OUT "trim_galore --paired --phred$phred --fastqc --illumina $sample1\.fastq.gz $sample2\.fastq.gz --output_dir ../fastq_trim\n";
	print OUT "bismark --bowtie2 --multicore 2 --phred$phred-quals --fastq -L 25 -N 1 $BismarkRefereDb -1 ../fastq_trim/$sample1\_val_1.fq.gz -2 ../fastq_trim/$sample2\_val_2.fq.gz -o ../bam\n";
	print OUT "filter_non_conversion --paired ../bam/$sample1\_val_1_bismark_bt2_pe.bam\n";
	print OUT "samtools sort ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam -o ../sortbam/$sample\_bismark_bt2_pe.sort.bam\n";
	print OUT "samtools index ../sortbam/$sample\_bismark_bt2_pe.sort.bam\n";
	print OUT "bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 10 --multicore 8 --paired-end --bedGraph --ignore 1 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam\n\n";
    }elsif(scalar(@read) == 1){	
  	my($sample,undef)=split /_1.fastq.gz|_R1.fastq.gz/,$read[0]; 	
    my($sample1,undef)=split /.fastq.gz/,$read[0];
    my $job_file_name = "$sample.pbs";
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");  
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $sample\n";
    print OUT "#PBS -q $queue\n";  # glean is free
    print OUT "#PBS -l nodes=$nodes:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime{$queue}\n";
    print OUT "#PBS -o $sample.log\n";
    print OUT "#PBS -e $sample.err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    
	print OUT "# fastq-dump --split-files --gzip $sample\n";
	print OUT "trim_galore --phred$phred --fastqc --illumina $sample1.fastq.gz --output_dir ../fastq_trim\n";
    print OUT "bismark --bowtie2 --phred$phred-quals --fastq -L 32 -N 1 --multicore 2 $BismarkRefereDb ../fastq_trim/$sample1\_trimmed.fq.gz -o ../bam\n";  
	print OUT "samtools sort ../bam/$sample\_se.bam -o ../sortbam/$sample.sort.bam\n";
	print OUT "samtools index ../sortbam/$sample.sort.bam\n";
	print OUT "bismark_methylation_extractor --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/$sample\_se.bam";
    }
   close(OUT);
}


sub USAGE{
print "------------------------------------------------------------------------------------\n\n
USAGE: perl $0 --input saminfo.txt --genome hg19  --server TSCC --submit\n";

print "saminfo.txt format: \n
SRR1648425_1.fastq.gz   SRR1648425_2.fastq.gz
SRR1648426_1.fastq.gz   SRR1648426_2.fastq.gz
SRR1648427_1.fastq.gz   SRR1648427_2.fastq.gz
-------------------------------------------------------------------\n";
}













