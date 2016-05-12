#!/usr/bin/perl -w
use strict;
use Cwd;
# Extract Pileup within specfic genomic regions (CDS, Exon)
# Set Pileup and GenomicInterval regions
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-03-16

chdir "C:\\Users\\shicheng\\Dropbox\\Project\\methylation\\monod\\analysis\\May52016-methHaplot";
my %sam;
my %hapinfo;
my %loc;
my %hap;
my %cpg;

die &USAGE if scalar @ARGV<4;
my $ct=shift @ARGV;
my $cp=shift @ARGV;
my $excl=shift @ARGV;
my $out=shift @ARGV;

# Test Hapinfo File
# my $ct="test1.txt";
# my $cp="test2.txt";
# my $out="rlt1.txt";


my ($ct_sam,undef)=split/\./,$ct;  # Should be keep consistent with Line 26: my ($sam,undef)=split/\./,$file;
my ($cp_sam,undef)=split/\./,$cp;  # Should be keep consistent with Line 26: my ($sam,undef)=split/\./,$file;

my @file=glob("*.hapInfo.txt");

foreach my $file(@file){
	my ($sam,undef)=split/\./,$file;
	open F,$file;
	while(<F>){
	chomp;
    next if /^\s+$/;
	my ($loc,$haptype,$count,$cpg)=split /\s+/;
	my $C_number = () = $haptype =~ /C/gi;		
		if($C_number>=2){                             # Criteron of HMH(High methylation Haplotype), Here at least two methylated C in the haplotype
		print "$sam\t$_\n";		
		$hapinfo{$loc}{$cpg}{$haptype}{$sam}=$count;
		$loc{$loc}=$loc;
		$cpg{$cpg}=$cpg;
		$hap{$haptype}=$haptype;
		$sam{$sam}=$sam;
		}
	}
}

open OUT,">$out.txt";
foreach my $loc(sort keys %hapinfo){
	my %Hapinfo;   # reformed hapinfo(hapsplit and hapmerge)
	foreach my $cpg(sort keys %{$hapinfo{$loc}}){
	foreach my $hap(sort keys %{$hapinfo{$loc}{$cpg}}){
		my @cpg=split/,/,$cpg;
		my $i=0;
		while($i<$#cpg){
			my $j=$i+1;
			while($j<=$#cpg-$i+1){
			my $string_hap=substr($hap,$i,$j);
			my $string_cpg=join(",",@cpg[$i..($i+$j-1)]);
			print "$string_hap\t$string_cpg";
			foreach my $sam(sort keys %sam){
			if($hapinfo{$loc}{$cpg}{$hap}{$sam}){
			$Hapinfo{$loc}{$cpg}{$hap}{$sam}+=$hapinfo{$loc}{$cpg}{$hap}{$sam};
			}
			}
			$j++;
			print "\n";
			}
			print "\n";
			$i++;				
		}	
	}
	}
	
	my @sam=sort keys %sam;
	my $header=join("\t",@sam);
	print OUT "LOC\tCpG\tHapType\t$header\n";
	foreach my $loc(sort keys %Hapinfo){
	foreach my $cpg(sort keys %{$Hapinfo{$loc}}){
	foreach my $hap(sort keys %{$Hapinfo{$loc}{$cpg}}){
		next if(! exists $Hapinfo{$loc}{$cpg}{$hap}{$ct_sam} || ! exists $Hapinfo{$loc}{$cpg}{$hap}{$cp_sam} || $Hapinfo{$loc}{$cpg}{$hap}{$ct_sam}<=1 || $Hapinfo{$loc}{$cpg}{$hap}{$cp_sam}<=1);  # Criteron of Minimum HMH in the Target Samples.
		print OUT "$loc\t$cpg\t$hap\t";	
		foreach my $sam(sort keys %sam){
		if(! exists  $Hapinfo{$loc}{$cpg}{$hap}{$sam}){
		print OUT "\t0";		
		}else{
		print OUT "\t$Hapinfo{$loc}{$cpg}{$hap}{$sam}";	
		}
		}
		print OUT "\n";
	}
	}
}
}


sub USAGE{
print "Usage: perl $0 M-Primary-Hapinfo M-Plasma-Hapinfo ExcludeDataBase-Hapinfo OutputFileName\n";
print "Identify the Potential Confounder of the Source Cancer Specific Methylation Haplotype\n";

print '
Format For: M-Primary-Hapinfo and M-Plasma-Hapinfo: 
chr1:100231312-100231328	CCT		1	100231316,100231324,100231328
chr1:100231312-100231328	CC		1	100231324,100231328
chr1:100231312-100231328	CCCC	4	100231312,100231316,100231324,100231328

chr1:100231312-100231328	CCC		1	100231316,100231324,100231328
chr1:100231312-100231328	CCCC	2	100231312,100231316,100231324,100231328
chr1:100231312-100231328	TTT		1	100231316,100231324,100231328
chr1:100231312-100231328	TT		1	100231324,100231328
chr1:100231312-100231328	T		2	100231312
chr1:100231312-100231328	CC		2	100231324,100231328
chr1:100231312-100231328	TTTT	4	100231312,100231316,100231324,100231328
chr1:100231312-100231328	C		1	100231328
';
print '
Format For: ExcludeDataBase-Hapinfo: Hapinfo File List to Excluded: 
STL001BL-01.hapInfo.txt
STL001FT-01.hapInfo.txt
STL001GA-01.hapInfo.txt
STL001LG-01.hapInfo.txt
STL001LV-01.hapInfo.txt
STL001PO-01.hapInfo.txt
STL001RV-01.hapInfo.txt
STL001SB-01.hapInfo.txt
STL001SG-01.hapInfo.txt
STL001SX-01.hapInfo.txt
';
}
