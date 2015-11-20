#!/usr/bin/perl
use strict;
my $dna=shift @ARGV;
my $rcdna= & reverse_complement_IUPAC($dna);
print "$rcdna\n";

sub reverse_complement_IUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}
