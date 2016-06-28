#!/usr/bin/perl
use strict;
my $file=shift @ARGV;

open F,$file;
chomp (my @bed=<F>);
close F;

foreach my $bed(@bed){
my @tmp=split/\s+/,$bed;
my $cor="$tmp[0]:$tmp[1]-$tmp[2]";
print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$cor\n";
}

