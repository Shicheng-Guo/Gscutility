#!/usr/bin/perl -w

# Merge Same Format Single Column files to One-Matrix
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-12-19

use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;
my @file=glob("R2*");

my %mhb;
my %pos;
my %sam;
foreach my $sam(@file){
$sam{$sam}=$sam;
open F,$sam;
while(<F>){
	chomp;
	my ($pos,$R2)=split /\t/;
	$mhb{$pos}{$sam}=$R2;
	$pos{$pos}=$pos;
	}	
}

######## Print Matrix-header
my @sam;
foreach my $sam(sort keys %sam){
	push(@sam,$sam);
}
my $head=join("\t",@sam);
print "\t$head\n";

foreach my $pos(sort keys %mhb){
 	print "$pos";
 	foreach my $sam (sort keys %sam){
 	if( ! exists $mhb{$pos}{$sam}||! defined $mhb{$pos}{$sam}||$mhb{$pos}{$sam}=~/NA/){
 	print "\tNA";			
 	}else{
 	my $R2=sprintf("%.3f",$mhb{$pos}{$sam});
 	print "\t$R2";
 	}	
 	}
 	print "\n";
}



