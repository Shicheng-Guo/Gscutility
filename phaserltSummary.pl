#!/usr/bin/perl

# Modify output file from phase to be easy to read and open with excel.
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-02-03

use strict;
use warnings;
die "perl phaserltSummary.pl\n" if scalar(@ARGV<1);
my $gencode_file=shift @ARGV;
my ($adjust,@files,$file,$dir,$n,@a);
@files=glob("*.out");
while(my $filename=shift @files){
	open FH, $filename;
	(my $GeneSymbol = $filename) =~ s/\.//ig; 
	while(<FH>){
		if (/BEGIN LIST_SUMMARY/){
			while (<FH>){
				my $lines=$_;
				my @temp=split /\s+/;
				for my $i(1..3){
					next if /END LIST_SUMMARY/;
					$a[$i]=$temp[2]  if $temp[1] eq $i;
				}
			}
		}
	}
}

print "$n haplotype were identfied\n";
