#!/usr/bin/perl

# Add the coordinate of introns to GENCODE GTF file based on transcript id.
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-02-03

use strict;
use warnings;
die "perl intronInsert2GTF GENCODE.gtf\n" if scalar(@ARGV<1);
my $gencode_file=shift @ARGV;

my @line;
open IN, $gencode_file or die "Can't open $gencode_file.\n";
@line=<IN>;
close IN;

my %gencode;
foreach my $i(0..($#line-1)){
    print "$line[$i]";
    my ($gene_id,$transcript_id);
    my %attribs = ();
    my ($chr1, $source1, $type1, $start1, $end1, $score1, $strand1, $phase1, $attributes1) = split("\t",$line[$i]);
    my $line2=$line[$i+1];
    my $j=$i+1;
    my ($chr2, $source2, $type2, $start2, $end2, $score2, $strand2, $phase2, $attributes2) = split("\t",$line[$j]);
    if(defined $type2){
    if($type2 eq 'exon' && $type1 eq 'exon'){
    my ($start3,$end3);
    if($strand1 eq '+'){
    $start3=$end1+1;
    $end3=$start2-1;
    }elsif($strand2 eq '-'){
    $start3=$end2+1;
    $end3=$start1-1;
    }
    $attributes1=~s/exon_number/intron_number/i;
    print "$chr1\t$source1\tintrion\t$start3\t$end3\t.\t$strand1\t$phase1\t$attributes1";
    }
}
}
