#!/usr/bin/perl
# add dbsnp153 to the file according to chr-pos-ref-alt in hg19 human genome
# contact: Shihcheng.Guo@gmail.com

use strict;
use Cwd;
my $dbsnp153="dbSNP153.hg19.vcf";
my %data;
#open F1,'<:gzip',$dbsnp153 || die "cannot open $dbsnp153";
open F1,$dbsnp153 || die "cannot open $dbsnp153";
while(<F1>){
next if /^#/;
my ($chr,$pos,$rs,$ref,$alt)=split/\s+/;
$ref=uc $ref;
$alt= uc $alt;
my $id="$chr-$pos-$ref-$alt";
$data{$id}=$rs;
}

my $input=shift @ARGV;
open F2, $input;
while(<F2>){
my @line=split/\s+/;
my $chr=$line[3];
my $pos=$line[4];
my $ref=uc $line[6];
my $alt=uc $line[5];
my $id="$chr-$pos-$ref-$alt";
if(defined $data{$id}){
$line[1]=$data{$id};
}
my $print= join("\t",@line);
print "$print\n";
}
