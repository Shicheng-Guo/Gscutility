#!/usr/bin/env perl

use strict;
use Cwd;
chdir getcwd;

my @input=glob("Chr*_450kMerge.txt");
foreach my $input(@input){
open F,$input;
my $output="$input.trans";
open OUT,">$output";
my (@rows, @transposed);
@rows = ();
@transposed = ();
while(<F>){
	# This is each row in your table
	chomp;
	my @row= split /\t/;
	push(@rows, \@row);
	}
	# Transfer #1
	for my $row (@rows){
		for my $column (0 .. $#{$row}) {
			push(@{$transposed[$column]}, $row->[$column]);
		}
	}
	# Transfer #2
	for my $new_row (@transposed) {
		for my $new_col (@{$new_row}) {
			print OUT $new_col,"\t";
		}
		print OUT "\n";
	}
}
