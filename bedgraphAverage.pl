#!/usr/bin/perl
use strict;
use Cwd;
use POSIX;
use List::Util qw(sum);
my $dir = getcwd;
open F,shift @ARGV;
my %wig;
while(<F>){
chomp;
my ($chr,$start,$end,$value)=split/\s+/;
push @{$wig{"$chr:$start-$end"}},$value;
}
foreach my $loci(sort keys %wig){
my ($chr,$start,$end)=split/:|-/,$loci;
my @line=@{$wig{$loci}};
print "$chr\t$start\t$end\tsum(@line)/($#line}+1)\n";
}
