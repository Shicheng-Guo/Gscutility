#!/usr/bin/perl -w

# Merge ChrM methlation level from bismark coverage output (Same Format Single Column files to One-Matrix)
# # 4: MF; #6: total coverage>10
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-12-19

use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;
my @file=glob("*chrM");

my %mhb;
my %pos;
my %sam;
foreach my $sam(@file){
$sam{$sam}=$sam;
open F,$sam;
while(<F>){
	chomp;
	my ($chr,$start,$end,$percent,$Mnum,$Tnum)=split /\t/;
        my $cpg="$chr:$start-$end";
	$pos{$cpg}=$cpg;
	next if $Tnum<10;
	$mhb{$cpg}{$sam}=$percent;
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



