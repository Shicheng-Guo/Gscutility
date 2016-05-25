#!/usr/bin/perl -w
use strict;
use Cwd;
# Extract Pileup within specfic genomic regions (CDS, Exon)
# Set Pileup and GenomicInterval regions
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-03-16

die USAGE() if scalar(@ARGV)<1;
# my $bamdir="/home/shg047/oasis/SALK/bam";
my $bamdir=shift @ARGV;
chdir $bamdir || die "$bamdir not avaibale\n";
my @file=glob("*bam");
# my $outdir="/home/shg047/oasis/Estellar2016/bam";
# STL011LI-01.chr16.sorted.clipped.bam.job;

my %sam;
foreach my $file(@file){
        my ($sam,$chr,undef)=split /\./,$file;
        $sam{$sam}++;
}

foreach my $sam(sort keys %sam){
	print "$sam\t$sam{$sam}\n";
}

sub USAGE{
        print "\nUSAGE: perl $0 bamDir\n\n";
        print "Example: perl $0 /media/LTS_60T/Dinh/BAM /media/LTS_60T/Shicheng/BAM\n";
}
