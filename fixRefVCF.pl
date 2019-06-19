#!/usr/bin/perl

use strict;
use warnings;

$ARGV[2] or die "use fixRefVCF.pl <GENOME_FASTA> <VCF> <FIX_VCF>\n";

my $genome = shift @ARGV;
my $orig_vcf = shift @ARGV;
my $new_vcf = shift @ARGV;

my %seq = ();
my ($chr, $pos, $ref, $obs);

warn "loading genome sequence from $genome\n";
open (my $fh, "<", $genome) or die "cannot read $genome\n";
while (<$fh>) {
    chomp;
    if (/>\w+(\d+)/) {
        $chr = $1;
    } else {
       $seq{$chr} .= $_;
    }
}
close $fh;

open (my $vh, "<", $orig_vcf) or die "cannot read $orig_vcf\n";
open (my $oh, ">", $new_vcf) or die "cannot write $new_vcf\n";
while (<$vh>) {
    if (/^#/) {
         print $oh $_;
    } else {
        chomp;
        my @var = split (/\t/, $_);
        $chr = $var[0];
        $pos = $var[1]; 
        $ref = $var[3];
        $obs = substr($seq{$chr}, $pos - 1, length $ref);
        if ($ref ne $obs) {
            warn "editing $chr:$pos $ref->$obs\n";
            $var[3] = $obs;
        }
        print $oh join "\t", @var;
        print $oh "\n";
    }
}
close $vh;
close $oh;
