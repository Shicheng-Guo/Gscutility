#!/usr/bin/perl

# cut chromosome from beginning to end with specified windows
# Usage: perl $0 chr21 1000
# Window=2000bp, Step=500bp

use strict;
use Cwd;
chdir getcwd;

my $CHR=shift @ARGV;
my $STEP=shift @ARGV;
my $chrLen="/oasis/tscc/scratch/shg047/db/hg19/hg19.chrom.sizes";
open F,$chrLen || die "cannot open $chrLen\n";
my %len;
while(<F>){
my($chr,$len)=split/\s+/;
$len{$chr}=$len;
}
my $start=1;
my $end=$len{$CHR};
while($start<($len{$CHR}-$STEP)){
$end=$start+2000;
my $id="$CHR:$start-$end";
print "$CHR\t$start\t$end\t$id\n";
$start=$start+$STEP+1;
}
