#!/usr/bin/perl
use strict;
use Cwd;
my $dir = getcwd;
chdir "/home/guosa/hpc/db/hg19/fa";

my $file=shift @ARGV;
my %chrosome;
my %cpgsnp;

open F,$file;
my $chr=split/\./,$file;
my $line;
while(<F>){
	next if />/;
	chomp;
	$line.=$_;
}
my @line=split //,$line;

my $bim="FinalRelease_QC_20140311_Team1_Marshfield.bim";

open F1,"$bim";
while(<F1>){
	chomp;
	my $cpgsnp=$_;
	my($chroms,$rs,undef,$pos,$ref,$alt)=split/\s+/;
	next if $chroms!=$chr;
	next if length($ref)>1;
	next if length($alt)>1;
    if("$ref$alt"=~/C/ && $line[$pos+1]=="G"){
      	$cpgsnp{$cpgsnp}= $cpgsnp;
    }elsif("$ref$alt"=~/G/ && $line[$pos-1]=="C"){
      	$cpgsnp{$cpgsnp}= $cpgsnp;
    }
}
foreach my $cpgsnp(sort keys %cpgsnp){
	print "$cpgsnp{$cpgsnp}\n";
}
print "Print CpG-SNP completed......\n";
