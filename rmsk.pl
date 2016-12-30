#!/usr/bin/perl -w

# reformat repeatmask format to bed format(rmsk.hg19 -> rmsk.hg19.bed)
# Good habit to creat parameter table for a project
# Contact: Shicheng Guo
# Version 1.3
# Update: 12/17/2016

use strict;
use Cwd;
use Getopt::Long;
chdir getcwd;

sub options{
	my $help;
	my $version;
	my $input;
	my $output;
    my $command_line=GetOptions (
    "help|h" => \$help,
    "version|v" => \$version,
    "input|i=s" => \$input,
    "ouput|o=s" => \$output
    );
    return($help,$version,$input,$output);
}	

my ($help,$version,$input,$output)=options();

open F,$input;
open OUT,">$output";
while(<F>){
next if /repClass/;
chomp;
my @line=split/\s/;
my $chr=$line[5];
my $start=$line[6];
my $end=$line[7];
my $id="$line[5]:$line[6]-$line[7]";
my $chain=$line[9];
my $type=$line[10];
my $group=$line[11];
print OUT "$chr\t$start\t$end\t$id\t$group\n";
}

