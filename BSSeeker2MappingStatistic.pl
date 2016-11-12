#!/usr/bin/perl -w

# BSSeeker2MappingStastic.pl: a perl script to collect the mapping information from BSSeeker2.
# Run the script to the Bam directory of the bismark
# Contact: Shicheng Guo
# Version 1.3

use strict;
use warnings;
use Cwd;

my $dir=getcwd;
chdir $dir;
print "\nenter to directory: $dir\nstarting reading mapping informations....\n";

# _bt2_PE_report.txt
# _bt2_SE_report.txt
# _bt2_PE_report.txt

my @file=glob("*bs_seeker2_log");

my $output="mapping.statistic.rlt.txt";
open OUT,">$output";
print OUT "Sample\tN(reads)\tN(base)\tN(multiplehit)\tN(Mappability)\tN(mappedbases\tN(non-map-reads)\tP(mCG)\tP(mCHG)\tP(mCHH)\n";
foreach my $file(@file){
open F,$file;
my ($sample,undef)=split /_trimmed|_val|.read1_val/,$file;
print OUT "$sample\t";
while(<F>){
chomp;
print OUT "$1\t" if /pairs:\s+(\d+)/i;
print OUT "$1\t" if /total:\s+(\w+)/;
print OUT "$1\t" if /hits:\s+(\w+)/;
print OUT "$1\t" if /post-filtering\):\s+(\d+.\d*%)/;
print OUT "$1\t" if /Mappability\s+=\s+(\d+.\d*%)/;
print OUT "$1\t" if /mapped reads:\s+(\d+)/;
print OUT "$1\t" if /Unmapped read pairs:\s+(\d+)/;
print OUT "$1\t" if /mCG\s+(\d+.\d*%)/;
print OUT "$1\t" if /mCHG\s+(\d+.\d*%)/;
print OUT "$1\t" if /mCHH\s+(\d+.\d*%)/;
}
print OUT "\n";
}

print "output file were saved in: $output\n\n";

