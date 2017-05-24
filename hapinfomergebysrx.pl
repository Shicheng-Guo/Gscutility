#!/usr/bin/perl -w
use strict;
use Cwd;
use Sort::Array qw/Sort_Table/;
chdir getcwd;
my $file=shift @ARGV;
my %SRA;
open F,$file;
while(<F>){
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
my %data;
print "reading $srr.hapInfo.txt...\n";
foreach my $srr(@{$SRA{$srx}}){
open F,"$srr.hapInfo.txt";
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

