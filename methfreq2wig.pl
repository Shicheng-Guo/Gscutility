#!/usr/bin/perl -w

# A perl script to collect zero-based coverage file to wig files so that it can be upload to UCSC.
# GSE52271, download methyfreq files and change them to bedgraph and bw
# Contact: Shihcheng.Guo@gmail.com
# Version 1.3
use strict;
use warnings;
use Cwd;

use strict;
my %tissue;
my @file=glob("GSM*.txt");

foreach my $file(@file){
open F2,$file;
my ($sample,undef)=split /\_/,$file;	
open OUT,">$sample..bedGraph";
my ($gsm,undef,$sam,undef)=split /\.|_/,$file;
print "track type=bedGraph name=\"$sam-$gsm\" visibility=full color=20,150,20 altColor=150,20,20 windowingFunction=mean\n";
while(<F2>){
chomp;
my @line=split /\t/;
next if ($line[2]+$line[3])<5;
my $chr=$line[0];
my $start=$line[1];
my $end=$line[1]+1;
my $mf=(100*$line[2])/($line[3]+$line[2]);
print OUT "$chr\t$start\t$end\t$mf\n";
}
system("sort -k1,1 -k2,2n $sample.WGBS.Estellar2016.bedGraph > $sample.WGBS.Estellar2016.sorted.bedGraph");
unlink "$sample.WGBS.Estellar2016.bedGraph";
system("bedGraphToBigWig $sample.WGBS.Estellar2016.sorted.bedGraph ~/oasis/db/hg19/hg19.chrom.sizes $sample.bw");
}

