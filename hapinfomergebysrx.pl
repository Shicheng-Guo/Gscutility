#!/usr/bin/perl -w
use strict;
use Cwd;
use Sort::Array qw/Sort_Table/;
use File::Copy;

chdir getcwd;
my $file=shift @ARGV;
my %SRA;
open F,$file;
while(<F>){
next if /sample_accession/;
next if /study_accession/;
chomp;
if(/(SRR\d+)/){
        my $SRR=$1;
        if(/(SRX\d+)/){
                my $SRX=$1;
                print "$SRX\t$SRR\n";
                push @{$SRA{$SRX}},$SRR;
                }
        }
}

my @SRX=sort keys %SRA;

foreach my $srx(@SRX){

if(scalar(@{$SRA{$srx}})<2){
my $srrname=shift @{$SRA{$srx}};
print "only one SRR file: $srrname point to $srx, direct copy will be conducted...\n";
copy("$srrname.hapInfo.txt","$srx.hapInfo.txt") or die "cannot open $srrname.hapInfo.txt or $srx.hapInfo.txt\n";
next;
}

my %data;
foreach my $srr(@{$SRA{$srx}}){
my($file,undef)=glob("$srr*hapInfo.txt");
print "reading $file\t...\n";
open F,"$file";
while(<F>){
	chomp;
	my ($mhb,$haptype,$number,$pos)=split /\t/;
        next if length($haptype)<1;
	my $key="$mhb\_$haptype\_$pos";
	$data{$key}+=$number
}
}
print "start writing $srx.hapInfo.txt ...\n";
open OUT,">$srx.hapInfo.txt";
foreach my $key(sort keys %data){
my ($mhb,$haptype,$pos)=split /_/,$key;
print OUT "$mhb\t$haptype\t$data{$key}\t$pos\n"; 
}
close OUT;
}

