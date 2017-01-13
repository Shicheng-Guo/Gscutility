#!/usr/bin/perl -w

# bismark to alignment single and pair-end fastq in same project
# Good habit to creat parameter table for a project
# Contact: Shicheng Guo
# Version 1.3
# Update: 1/17/2017

use strict;
use Cwd;
use Getopt::Long;

my ($input,$submit,$genome,$server,$help,$queue,%ppn,%multicore,%walltime);

## default parameters
$input ="saminfo.txt";
$submit='nonsubmit';
$genome="hg19";
$server="TSCC";

## parse parameters
GetOptions ( "input=s"   => \$input,          # string
             "genome=s" => \$genome,          # string   
             "server=s" => \$server,          # flag
             "help" => \$help,                # help
	     "queue=s" => \$queue,              # queue(hotel,pdafm,condo)
			 )          
or die("Error in command line arguments\n");

print @ARGV;

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
    
    %walltime=(
    "hotel" => "168:00:00",
    "condo" => "8:00:00",
    "pdafm" => "72:00:00",
    );
    
    %ppn=(
    hotel => "16",
    pdafm => "32",
    glean => "16",
    condo => "16",
    );

    %multicore=(
    hotel => "6",
    pdafm => "12",
    glean => "6",
    condo => "6",
    );
	
	my $nodes=1;
    my $curr_dir = $dir;
    
    
    if(scalar(@read) eq 2){
    my($sample,undef)=split /_1.fastq.gz|_R1.fastq.gz/,$read[0]; 	
    my($sample1,undef)=split /.fastq.gz/,$read[0];
    my($sample2,undef)=split /.fastq.gz/,$read[1];
    my $job_file_name = "$sample.pbs";
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");   
    
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
    print OUT "# fastq-dump --split-files --gzip $sample\n";
	print OUT "trim_galore --paired --phred$phred --fastqc --illumina $sample1\.fastq.gz $sample2\.fastq.gz --output_dir ../fastq_trim\n";
	print OUT "bismark --bowtie2 --multicore $multicore{$queue} --phred$phred-quals --fastq -L 25 -N 1 $BismarkRefereDb -1 ../fastq_trim/$sample1\_val_1.fq.gz -2 ../fastq_trim/$sample2\_val_2.fq.gz -o ../bam\n";
	print OUT "filter_non_conversion --paired ../bam/$sample1\_val_1_bismark_bt2_pe.bam\n";
	print OUT "rm ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_removed_seqs.bam\n";
	print OUT "samtools sort -@ 8  ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam -o ../sortbam/$sample\_bismark_bt2_pe.sort.bam\n";
	print OUT "samtools index ../sortbam/$sample\_bismark_bt2_pe.sort.bam\n";
	print OUT "bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 1 --multicore 8 --paired-end --bedGraph --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../bam/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bam\n";
    print OUT "coverage2cytosine --merge_CpG ../methyfreq/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov.gz --genome_folder $BismarkRefereDb -o ../methyfreq/$sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bismark.CpG_merge.cov";
	}elsif(scalar(@read) == 1){	
  	my($sample,undef)=split /_1.fastq.gz|_R1.fastq.gz/,$read[0]; 	
    my($sample1,undef)=split /.fastq.gz/,$read[0];
    my $job_file_name = "$sample.pbs";
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");  
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
	print OUT "# fastq-dump --split-files --gzip $sample\n";
	print OUT "trim_galore --phred$phred --fastqc --illumina $sample1.fastq.gz --output_dir ../fastq_trim\n";
    print OUT "bismark --bowtie2 --phred$phred-quals --fastq -L 32 -N 1 --multicore $multicore{$queue} $BismarkRefereDb ../fastq_trim/$sample1\_trimmed.fq.gz -o ../bam\n";  
	print OUT "filter_non_conversion --paired ../bam/$sample1\_val_1_bismark_bt2_pe.bam\n";
	print OUT "rm ../bam/$sample1\_val_1_bismark_bt2_.nonCG_removed_seqs.bam\n";
	print OUT "samtools sort ../bam/$sample\_se.bam -o ../sortbam/$sample.sort.bam\n";
	print OUT "samtools index ../sortbam/$sample.sort.bam\n";
	print OUT "bismark_methylation_extractor --single-end --bedGraph --cutoff 1 --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../bam/$sample\_se.bam\n";
    print OUT "coverage2cytosine --merge_CpG ../methyfreq/$sample\_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov.gz --genome_folder $BismarkRefereDb -o ../methyfreq/$sample\_val_1_bismark_bt2_pe.nonCG_filtered.bismark.CpG_merge.cov";
	}
   close(OUT);
}
