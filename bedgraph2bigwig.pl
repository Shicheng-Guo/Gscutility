#!/usr/bin/perl
use strict;

my $file=shift @ARGV;

my $rmheadSort="(awk 'NR!=1' $file| sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4}') > $file.sort.tmp";
my $bedgraph2bigwig="bedGraphToBigWig 6-T-4.sorted.clipped_CpG.sort.bedGraph ~/oasis/db/hg19.chrom.sizes 6-T-4.sorted.clipped_CpG.sort.bedGraph.bw";
system($rmheadSort)
system($bedgraph2bigwig)
