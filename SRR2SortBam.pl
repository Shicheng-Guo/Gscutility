#!/usr/bin/perl -w

# A perl script to build config file for SRS Bam Merge
# Contact: Shihcheng.Guo@Gmail.com
# Version 1.3
# Go to http://sra.dnanexus.com/studies/SRP028600/samples
# Select SRS and Click Related RUNS then get the Table as the input

use strict;
use warnings;
use Cwd;
my $curr_dir = getcwd;
my $file=shift @ARGV;
my $submit=shift @ARGV;
my %SRA;

my @SRRBAM=glob("*bam");

open F,$file;
while(<F>){
chomp;
if(/(SRR\d+)/){
        my $SRR=$1;
        if(/(SRX\d+)/){
                my $SRX=$1;
                my @SrrBam=grep(/$SRR/,@SRRBAM);
				next if scalar(@SrrBam)<1;
                print "$SrrBam[0]\n";
                push @{$SRA{$SRX}},$SrrBam[0];
                }
        }
}

mkdir "../mergeBam" if ! -e "../mergeBam";
mkdir "../methyfreq" if ! -e "../methyfreq";
mkdir "../SortMergeBam" if ! -e "../SortMergeBam";


foreach my $SRX(sort keys %SRA){
	
    my $id=$SRX;
    my $job_file_name = $id . ".job";
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -q glean\n";
    print OUT "#PBS -l nodes=1:ppn=2\n";
    print OUT "#PBS -l walltime=168:00:00\n";
    print OUT "#PBS -o ".$id.".log\n";
    print OUT "#PBS -e ".$id.".err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shicheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $curr_dir\n";
    my $srrList=join(" ",@{$SRA{$SRX}});
    print OUT "samtools merge -@ 2 ../mergeBam/$SRX.bam $srrList\n";
    print OUT "samtools sort ../mergeBam/$SRX.bam ../SortMergeBam/$SRX.sortn.bam\n";
    print OUT "samtools index ../SortMergeBam/$SRX.sortn.bam \n";
	#   my $is_pairedEnd= & is_bamPE("@{$SRA{$SRX}}[0].bam");
    my $is_pairedEnd = "--paired-end";
    print OUT "bismark_methylation_extractor $is_pairedEnd --bedGraph --multicore 2 --ignore 3 --ignore_3prime 3 --ignore_r2 5 --ignore_3prime_r2 5 --gzip --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../SortMergeBam/$SRX.sortn.bam\n";
    print OUT "samtools sort ../mergeBam/$SRX.bam ../SortMergeBam/$SRX.sortc.bam\n";
    print OUT "samtools index ../SortMergeBam/$SRX.sortc.bam \n";

    if($submit eq 'submit'){
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
close OUT;




