#!/usr/bin/perl -w

# estimate the length of the Introns/intron/transcript from GENCODE GTF file
# remember to run intronInsert.pl first and then run the script to the output of introninsert.pl
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-02-8
use strict;
use warnings;
die "perl intronEstimateFromGTF.pl GENCODE.gtf\n" if scalar(@ARGV<1);
my $gencode_file=shift @ARGV;

open IN, $gencode_file or die "Can't open $gencode_file.\n";
my %gencode;
while(<IN>){
        next if ! /transcript|exon|intron/;
        my ($gene_id,$transcript_id);
        my %attribs = ();
        my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t");
    my @add_attributes = split(";", $attributes);
    my (@c_type,@c_value);
    foreach my $attr (@add_attributes) {
        next unless $attr =~ /^\s*(.+)\s(.+)$/;
        my $c_type  = $1;
        my $c_value = $2;

        push(@c_type,$c_type);
        push(@c_value,$c_value);
    }
    my $len=$end-$start;
    push @{$gencode{$c_value[0]}{$c_value[1]}{$type}},$len;
}

        print "\tgene_id\ttransLen\ttransNum\texonLen\texonNum\tintronLen\tintronNum\n";

foreach my $gene_id(sort keys %gencode){

        my $transLen=0;
        my $exonLen=0;
        my $transNum=0;
        my $intronLen=0;
        my $exonNum=0;
        my $intronNum=0;

        foreach my $transcript_id(sort keys %{$gencode{$gene_id}}){
        $transNum++;
        foreach my $trans_len(@{$gencode{$gene_id}{$transcript_id}{'transcript'}}){
        $transLen+=$trans_len;
        }
        foreach my $exon_len(@{$gencode{$gene_id}{$transcript_id}{'exon'}}){
        $exonLen+=$exon_len;
        $exonNum++;
        }
        foreach my $intron_len(@{$gencode{$gene_id}{$transcript_id}{'intron'}}){
        $intronLen+=$intron_len;
        $intronNum++;
        }
     }
   $gene_id=~s/"//g;
   print "$gene_id\t$transLen\t$transNum\t$exonLen\t$exonNum\t$intronLen\t$intronNum\n";
}

