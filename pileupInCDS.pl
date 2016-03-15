#!/usr/bin/perl -w
# Extract Pileup within specfic genomic regions (CDS, Exon)
# Set Pileup and GenomicInterval regions
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-03-15

use strict;
die &USAGE if scalar(@ARGV)<1;
my $input=shift @ARGV;
my $targetFile=shift @ARGV;
my $flankingLen = 0;
my $binSize= 10000;
my %exon_bin_table;
# my $targetFile = "/home/zhl002/genomeDB/mm9/CCDS.20120710.UCSC.chr.bed";
open(INFILE, "$targetFile") || die("Error in opening targetFile!\n");
while(my $line = <INFILE>){
        chop($line);
        my ($chr,$start,$end) = split(/\t/, $line);
        $chr =~ s/chr//;
        my $index = int(($start+$end)/2/$binSize);
        push(@{$exon_bin_table{$chr}->{$index-1}}, "$chr:$start:$end");
        push(@{$exon_bin_table{$chr}->{$index}}, "$chr:$start:$end");
        push(@{$exon_bin_table{$chr}->{$index+1}}, "$chr:$start:$end");
}
close(INFILE);

open F,$input;
while(my $line = <F>){
        my @fields = split(/\t/, $line);
        my $chr = $fields[0];
        $chr =~ s/chr//;
        $chr =~ s/\.//;
        $chr =~ s/fa//;
        my $index = int($fields[1]/$binSize);
        next if(!$exon_bin_table{$chr}->{$index});
        my @candidate_exons = @{$exon_bin_table{$chr}->{$index}};
        if($exon_bin_table{$chr}->{$index-1})
        {
                push(@candidate_exons,@{$exon_bin_table{$chr}->{$index-1}});
        }
        if($exon_bin_table{$chr}->{$index+1})
        {
                push(@candidate_exons,@{$exon_bin_table{$chr}->{$index+1}});
        }
        if($exon_bin_table{$chr}->{$index-2})
        {
                push(@candidate_exons,@{$exon_bin_table{$chr}->{$index-2}});
        }
        if($exon_bin_table{$chr}->{$index+2})
        {
                push(@candidate_exons,@{$exon_bin_table{$chr}->{$index+2}});
        }
        foreach my $exon_info (@candidate_exons){
                my ($chr, $start, $end) = split(/:/, $exon_info);
                if($fields[1] >= $start-$flankingLen && $fields[1] <=$end+$flankingLen){
                        print $line;
                        last;
                }
        }
}


sub USAGE{
print "\nUSAGE: $0 Pileup_File CDS_UCSC_File > Output.txt\n";
print "\nExtract Pileup within specfic genomic regions (CDS, Exon)\nAssign Pileup and GenomicInterval regions files as the paramters\n";
print '
CCDS File Format:
chr1    134212806       134213049       CCDS15283.1_exon_0_0_1_134212807_f      0       +
chr1    134221529       134221650       CCDS15283.1_exon_1_0_1_134221530_f      0       +
chr1    134224273       134224425       CCDS15283.1_exon_2_0_1_134224274_f      0       +
chr1    134224707       134224773       CCDS15283.1_exon_3_0_1_134224708_f      0       +
';	
print '
Pileup File Format:
chrM    1       G       G       255     0       56      309     ^].^].^D.^].^].^].^].^].^].^].^].^?,^A,^],^Q,^X,^],^Q,^],^],^],^],^],^],^],^],^M,^],^],^],^],^A,^],^],^],^],^],^],^],^],^V,^I,^],^],^],^],^K,^U,^Q,^U,^],^],^],^],^],^Q,^],^]
chrM    2       T       T       255     0       56      310     ...........,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
chrM    3       T       T       255     0       56      312     ...........,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
chrM    4       A       A       255     0       56      315     ...........,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
chrM    5       A       A       255     0       56      316     ...........,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
chrM    6       T       T       255     0       56      318     ....C......,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
';
}
