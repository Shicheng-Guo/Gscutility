#!/usr/bin/perl -w

# Transfer Bam to Fastq with samtools command
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-02-19
use strict;
use Cwd;
chdir getcwd;

# my $region=shift @ARGV;
& printUsage if scalar(@ARGV)<2;
my $region=shift @ARGV;       # chr5:112,039,140-112,047,142
my $r2File = shift @ARGV;  # R2 Matrix = "APC.chr5.rsq"
my $cpg_position_file=shift @ARGV; # "/home/kunzhang/HsGenome/hg19/HsGenome19.CpG.positions.txt";   # genome-miner
my @region=locationSplit($region);
open F,$cpg_position_file;
chomp(my @loci=<F>);
close F;
my $Loci=\@loci;
my @lociByRegion=lociPickUpByRegion($region,$Loci);
my %R2;
foreach my $loc1(@lociByRegion){
foreach my $loc2(@lociByRegion){
$R2{$loc1}{$loc2}="NA";
}
}

open F,$r2File;
my @R2loci;
while(<F>){
chomp;
if (/^\s/){
@R2loci=split /\t/;
next;
}
my @line=split /\s+/;
foreach my $tmp(1..$#line){
$R2{$line[0]}{$R2loci[$tmp]}=$line[$tmp];
}
}

foreach my $loc1(@lociByRegion){
my @tmp;
foreach my $loc2(@lociByRegion){
push(@tmp,$R2{$loc1}{$loc2});
}
my $print=join("\t",@tmp);
print "$print\n";
}

sub locationSplit{
my $location=shift;
$location=~s/,//g;
my ($chr, $start, $end) = split /[:-]/, $location;
return($chr,$start,$end);
}

sub lociPickUpByRegion{
	my $region=shift;
	my $Loci=shift;
	my @loci=@$Loci;
	my @lociPickup;
	my($chr1,$start1,$end1)=locationSplit($region);
	foreach my $loci(@loci){
	my ($chr2,$start2)=split/\t/,$loci;
	if ($chr2 eq $chr1 and $start2>=$start1 and $start2<= $end1){
	push(@lociPickup,$start2);
	}
	}	
	# @lociPickup=sort(@lociPickup);
	return(@lociPickup)
}

sub printUsage{
        print " Usage:\n";
		print " perl $0 <chr:start-end> <R-square matrix file> <human_cpg_position_file>\n";
		print " For example:\n perl $0 chr10:10000873-10001472 APC.chr5.rsq /home/kunzhang/HsGenome/hg19/HsGenome19.CpG.positions.txt > result.txt\n\n";
		print "-----------------------------------------------------------------------------\n";
		print " APC.chr5.rsq format:\n";
		print "           10000873  10001472 10001596\n";
		print "10000873    0.3         0.8      0.1   \n";
		print "10001472    0.5         0.4      0.9   \n";
		print "10001596    0.2         0.1      0.8   \n";
    	print "-----------------------------------------------------------------------------\n";
		print " human_cpg_position_file format:\n";
		print "chr1     10000873\n";
		print "chr1     10001472\n";
		print "chr1     10001596\n";
    	print "-----------------------------------------------------------------------------\n";
        exit 0;
}




