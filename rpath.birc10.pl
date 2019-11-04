#!/usr/bin/perl
use strict;
use Cwd;
my $dir=getcwd;
$dir=~s/mnt/\/mcrfnas2/ig;
$dir=~s/\//\//g;
print "$dir\n";
