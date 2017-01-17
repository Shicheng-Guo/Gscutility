#!/usr/bin/perl -w

# Separate/split refGene to Enhancer, Promoter, UTR, Exon, Intron subsection bed file. 
# Contact: Shicheng Guo
# Version 1.3
# Update: 2017-02-08
# Input: refGene.txt download from UCSC

use strict;
use warnings;
die "perl intronEstimateFromGTF.pl refGene.txt > refGene.bed\n" if scalar(@ARGV<1);
my $refGene=shift @ARGV;

open IN, $refGene or die "Can't open $refGene\n";
my %refGene;
while(<IN>){
	next /"random"|"hap"|"chrUN"/;
    my ($bin,$NM,$chr,$strand,$txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds,undef,$genesymbol,undef) = split /\t/;
    my @exonStarts = split(",", $exonStarts);
    my @exonEnds = split(",", $exonEnds);
    if($strand eq "+"){
    # gene in positive strand	
    my @enhancer=($txStart-4000,$txStart-2001);
    my @promter=($txStart-2000,$txStart-1);
    print "$chr\t$enhancer[0]\t$enhancer[0]\t$strand\t$NM\t$genesymbol\tEnhancer\n";
    print "$chr\t$promter[0]\t$promter[0]\t$strand\t$NM\t$genesymbol\tPromoter\n";
    my $exonrank;
    foreach my $i(0..($#exonStarts-1)){
    $exonrank++;
    print "$chr\t$exonStarts[$i]\t$exonEnds[$i]\t$strand\t$NM\t$genesymbol\tExon$exonrank\n";
    print "$chr\t$exonEnds[$i]\t$exonStarts[$i+1]\t$strand\t$NM\t$genesymbol\tIntron$exonrank\n";
    }
    $exonrank++;
    print "$chr\t$exonStarts[$#exonStarts]\t$exonEnds[$#exonStarts]\t+\t$NM\t$genesymbol\tExon$exonrank\n";
    }else{
    # gene in negative strand
    my @enhancer=($txEnd+2001,$txEnd+4000);
    my @promter=($txEnd-1,$txEnd-2000);
    print "$chr\t$enhancer[0]\t$enhancer[0]\t$strand\t$NM\t$genesymbol\tEnhancer\n";
    print "$chr\t$promter[0]\t$promter[0]\t$strand\t$NM\t$genesymbol\tPromoter\n";
    my $exonrank;
    foreach my $i(($#exonStarts)..1){
    $exonrank++;
    print "$chr\t$exonEnds[$i]\t$exonStarts[$i]\t+\t$NM\t$genesymbol\tExon$exonrank\n";
    print "$chr\t$exonStarts[$i]\t$exonEnds[$i-1]\t+\t$NM\t$genesymbol\tIntron$exonrank\n";
    }
    $exonrank++;
    print "$chr\t$exonEnds[0]\t$exonStarts[0]\t+\t$NM\t$genesymbol\tExon$exonrank\n";
    }        
}

        
