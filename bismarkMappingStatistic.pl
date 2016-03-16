#!/usr/bin/perl -w

# BismarkAlignment.pl: a perl script to collect the mapping information from bismark.
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
my @file=glob("*bt2*report.txt");

my $output="mapping.statistic.rlt.txt";
open OUT,">$output";
print OUT "Sample\tN(reads)\tN(mapped)\tP(mapping)\tN(C)\tN(MCPG)\tN(MCHG)\tN(MCHH)\tN(UCPG)\tN(UCHG)\tN(UCHH)\tP(MCPG)\tP(MCHG)\tP(MCHH)\n";
foreach my $file(@file){
open F,$file;
my ($sample,undef)=split /_trimmed|_val|.read1_val/,$file;
print OUT "$sample\t";
while(<F>){
chomp;
print OUT "$1\t" if /Sequence.*analysed in total:\s+(\d+)/i;
print OUT "$1\t" if /unique best hit:\s+(\w+)/;
print OUT "$1\t" if /from the different alignments:\s+(\w+)/;
print OUT "$1\t" if /Mapping efficiency:\s+(\d+.\d*%)/;
print OUT "$1\t" if /Total number of C\'s analysed:\s+(\w+)/;
print OUT "$1\t" if /Total methylated C\'s in CpG context:\s+(\d+)/;
print OUT "$1\t" if /Total methylated C\'s in CHG context:\s+(\d+)/;
print OUT "$1\t" if /Total methylated C\'s in CHH context:\s+(\d+)/;
print OUT "$1\t" if /Total unmethylated C\'s in CpG context:\s+(\d+)/;
print OUT "$1\t" if /Total unmethylated C\'s in CHG context:\s+(\d+)/;
print OUT "$1\t" if /Total unmethylated C\'s in CHH context:\s+(\d+)/;
print OUT "$1\t" if /C methylated in CpG context:\s+(\d+.\d*%)/;
print OUT "$1\t" if /C methylated in CHG context:\s+(\d+.\d*%)/;
print OUT "$1\t" if /C methylated in CHH context:\s+(\d+.\d*%)/;
}
print OUT "\n";
}
