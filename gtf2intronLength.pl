#!/usr/bin/perl -w

# estimate the length of the Introns from GENCODE GTF file
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-01-19
# Contact: Shicheng Guo<scguo@ucsd.edu>; Kun Zhang<kzhang@eng.ucsd.edu>
use strict;
use warnings;
die "perl intronEstimateFromGTF.pl GENCODE.gtf\n" if scalar(@ARGV<1);
my $gencode_file=shift @ARGV;

open IN, $gencode_file or die "Can't open $gencode_file.\n";
my %gencode;
while(<IN>){
        next if ! /transcript|exon/;
        my ($gene_id,$transcript_id);
        my %attribs = ();
        my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t");
    my @add_attributes = split(";", $attributes);
    my (@c_type,@c_value);
    foreach my $attr (@add_attributes ) {
        next unless $attr =~ /^\s*(.+)\s(.+)$/;
        my $c_type  = $1;
        my $c_value = $2;
        push(@c_type,$c_type);
        push(@c_value,$c_value);
    }
    my $len=$end-$start;
    push @{$gencode{$c_value[0]}{$c_value[2]}{$type}},$len;
}

foreach my $gene_id(sort keys %gencode){
        my ($transLen,$exonLen,$transNum);
        foreach my $transcript_id(sort keys %{$gencode{$gene_id}}){
                $transNum++;
                foreach my $trans_len(@{$gencode{$gene_id}{$transcript_id}{'transcript'}}){
                $transLen+=$trans_len;
                }
                foreach my $exon_len(@{$gencode{$gene_id}{$transcript_id}{'exon'}}){
                $exonLen+=$exon_len;

        }
     }
   my $intronLen=$transLen-$exonLen;
   print "$gene_id\t$transLen\t$exonLen\t$intronLen\t$transNum\n";
}
