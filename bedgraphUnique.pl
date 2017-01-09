#!/usr/bin/perl -w

# A foundmental script to merge redundant bedgraph to unique bedgraph
# Contact: Shihcheng.Guo@Gmail.com
# Version 1.3
# Go to http://sra.dnanexus.com/studies/SRP028600/samples
# input: redundant bedgraph
# output: non-redundant bedgraph
# Run the script in the fold of coverage files created by bismark alignmetor.

use strict;
use warnings;
use Cwd;
use Getopt::Long;

my $input=shift @ARGV;

open F,$input;
my %out;
while(<F>){
my ($chr,$start,$end,$mf,$me,$ume)=split/\s+/;
my $pos="$chr:$start-$end";
$out{$pos}{'me'}+=$me;
$out{$pos}{'ume'}+=$ume;
}

foreach my $pos(sort keys %out){
my $me=$out{$pos}{'me'};
my $ume=$out{$pos}{'ume'};
if(($me+$ume)>=1){
my $mf=sprintf("%.3f",($me/($me+$ume+0.01)));

my ($chr,$start,$end)=split/[:-]/,$pos;
print "$chr\t$start\t$end\t$mf\n";
}
}
