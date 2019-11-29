#!/usr/bin/perl
# add dbsnp153 to the file according to chr-pos-ref-alt in hg19 human genome
# contact: Shihcheng.Guo@gmail.com

use strict;
use Cwd;
my $dbsnp153="/home/guosa/hpc/db/dbSNP153/dbSNP153.hg19.vcf.gz";
my %data;
open F1,'<:gzip',$dbsnp153;
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
my @line=split/\s/+
my $chr=$line[3];
my $pos=$line[4];
my $ref=uc $line[5];
my $alt=uc $line[6];
my $id="$chr-$pos-$ref-$alt";
if(defined $data{id}){
$line[1]=$data{$id}
}
my $print= join("\t",@line);
}
