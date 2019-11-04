#!/usr/bin/perl
use strict;
use Cwd;
my $dir=getcwd;

$dir=~s/mnt/\\mcrfnas2/ig;
$dir=~s/\//\\/g;

print "$dir\n";

/gpfs/home/guosa/hpc/rheumatology/RA/he2020/impute/R3

