#!/usr/bin/perl -w

# Merge methlation level from bismark coverage output (Same Format Single Column files to One-Matrix)
# # 4: MF; #6: total coverage>10
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-12-19

# Following R script to make boxplot to check the variation
# mitometh<-read.table("C:\\Users\\shicheng\\Downloads\\CHR.Methylation.txt",head=T,row.names=1,as.is=T)
# colnames(mitometh)<-unlist(lapply(colnames(mitometh),function(x) substr(x,1,6)))
# rownames(mitometh)<-lapply(rownames(mitometh),function(x) unlist(strsplit(x,"[:-]"))[2])
# boxplot(t(mitometh),outline = F,horizontal = T)
# par(las=2,cex=0.6,mar=c(6,6,1,4))
# boxplot(t(mitometh),outline = T,horizontal = T)
# par(las=2,cex=0.6,mar=c(6,6,1,4))
# boxplot((mitometh),outline = T,horizontal = T)

use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;
my @file=glob("*bismark.cov");

my %mf;
my %pos;
my %sam;
foreach my $sam(@file){
my ($sam)=split/\_R1/,$sam;
$sam{$sam}=$sam;
open F,$sam;
while(<F>){
	chomp;
	my ($chr,$start,$end,$mf)=split /\t/;
  my $cpg="$chr:$start";
	$pos{$cpg}=$cpg;
	$mf{$cpg}{$sam}=$mf;
	}	
}

######## Print Matrix-header
my @sam;
foreach my $sam(sort keys %sam){
	push(@sam,$sam);
}
my $head=join("\t",@sam);
print "\t$head\n";

foreach my $pos(sort keys %mf){
 	print "$pos";
 	foreach my $sam (sort keys %sam){
 	if( ! exists $mf{$pos}{$sam}||! defined $mf{$pos}{$sam}||$mf{$pos}{$sam}=~/NA/){
 	print "\tNA";			
 	}else{
 	my $R2=sprintf("%.3f",$mf{$pos}{$sam});
 	print "\t$R2";
 	}	
 	}
 	print "\n";
}
