#!/usr/bin/perl -w

# estimate the length of the Introns from GENCODE GTF file
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-01-28
# Contact: Shicheng Guo<scguo@ucsd.edu>; Kun Zhang<kzhang@eng.ucsd.edu>

use strict;
use warnings;
die "perl ~/bin/table2wiki.pl input.txt\n" if scalar(@ARGV<1);
my $input=shift @ARGV;

open F,$input;
print "{| class=\"wikitable\" style=\"text-align: right; color: red;font-size:75%\"\n";
while(<F>){
chomp;
my @line=split /\t/;
my $tmp=join("||",@line);
print "| $tmp\n|-\n";
}
print"|}\n";
