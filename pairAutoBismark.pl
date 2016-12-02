#!/usr/bin/perl -w

# BismarkAlignment.pl: a perl script to run bismark for pair-end fastq.
# Run the script trim_galore, adaptor remove, bismark alignment and methylfreq and so on
# Contact: Shicheng Guo
# Version 1.3

use strict;
use warnings;
use Cwd;

my $dir=getcwd;
chdir $dir;

my $file=shift @ARGV;
my $phred=shift @ARGV;
my $submit=shift @ARGV;

my %SRA;
open F,$file;
while(<F>){
chomp;
if(/(SRR\d+)/){
	my $SRR=$1;
	if(/(SRX\d+)/){
		my $SRX=$1;
		print "$SRR\t$SRX\t$SRR.bismark.job\n";
		push @{$SRA{$SRX}},$SRR;
		}
	}
}


mkdir "../fastq_trim" if ! -e "../fastq_trim";
mkdir "../bam" if ! -e "../bam";
mkdir "../bedgraph" if ! -e "../bedgraph";
mkdir "../sortbam" if ! -e "../sortbam";
mkdir "../hapinfo" if ! -e "../sortbam";

foreach my $SRX(sort keys %SRA){
	foreach my $SRR (@{$SRA{$SRX}}){
        my $sample="$SRR";
	my $sample1="$SRR\_1";
	my $sample2="$SRR\_2";
	open OUT,">$SRR.bismark.job";
	# chomp(my $phredcheck=`perl /home/shg047/bin/checkphred.pl $sample1.fastq.gz`);
	# my ($phred)=split /\s+/,$phredcheck;

	print OUT "#!/bin/csh\n";
	print OUT "#PBS -q hotel\n";
	print OUT "#PBS -l nodes=1:ppn=8\n";
	print OUT "#PBS -l walltime=168:00:00\n";
	print OUT "#PBS -o $SRR\.log\n";
	print OUT "#PBS -e $SRR\.err\n";
	print OUT "#PBS -V\n";
	print OUT "#PBS -M shihcheng.guo\@gmail.com \n";
	print OUT "#PBS -m abe\n";
	print OUT "#PBS -A k4zhang-group\n";
	print OUT "cd $dir\n";
	# print OUT "fastq-dump --split-files --gzip $SRR\n";
	print OUT "trim_galore --paired --phred$phred --fastqc --illumina $sample1.fastq.gz $sample2.fastq.gz --output_dir ../fastq_trim\n";
	print OUT "bismark --bowtie2 --phred$phred-quals --multicore 2 --fastq -L 32 -N 1 /home/shg047/db/hg19/bismark/ -1 ../fastq_trim/$sample1\_trimmed.fq.gz -2 ../fastq_trim/$sample2\_trimmed.fq.gz -o ../bam\n";
	print OUT "samtools sort ../bam/$sample1\_trimmed.fq_bismark_bt2.bam -o ../sortbam/$sample\_trimmed.fq_bismark_bt2.sort.bam\n";
	print OUT "samtools index ../sortbam/$sample\_trimmed.fq_bismark_bt2.sort.bam\n";
	print OUT "bismark_methylation_extractor --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/$sample1\_trimmed.fq_bismark_bt2.bam";
	close OUT;
	if ($submit eq "submit"){
	system("qsub $SRR.bismark.job");
	}		
	}
}


sub print_helpfile{
  print << "HOW_TO";
  This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.

DESCRIPTION
USAGE: bismark [options] <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}
HOW_TO
}












