#!/usr/bin/perl -w

# Convert bed3 to bed4
# Contact: Shicheng Guo
# Version 1.3
# Update: Jan/19/2016

use strict;
use Cwd;
my $dir=getcwd;

my $bed3=shift @ARGV;
my $extension=200;
open F,$bed3;
while(<F>){
chomp;
my ($chr,$start,$end)=split/\t/;	
my $start2=$start-$extension;
my $end2=$end+$extension;
my $cor="$chr:$start2-$end2";
print "$chr\t$start2\t$end2\t$cor\n";
}
