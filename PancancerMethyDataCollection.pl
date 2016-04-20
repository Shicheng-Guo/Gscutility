#!/usr/bin/perl -w
#version 3.1
#update: mutiple cancer + latest version + special cpgsite 
#update: filerecords
#updata: April 20 2016
# jhu-usc.edu_LUSC.HumanMethylation27.3.lvl-3.TCGA-BG-A0MS-01A-11D-A10Q-05.txt 
# jhu-usc.edu_LGG.HumanMethylation450.11.lvl-3.TCGA-HT-A617-01A-11D-A29T-05.txt

use strict;
use Cwd;
chdir getcwd;

# my $annoFile=$ARGV[1]; # 450K annotation file
my $annoFile="/home/shg047/oasis/db/GPL13534.sort.bed";
# open OUT1, ">PancancerFileList.txt";
# open OUT2, ">PancancerGeneList.txt";
# open OUT3, ">PancancerCancerList.txt";

my (%PancancerSampleList,%PancancerGeneList,%PancancerCancerList,%MethDataMatrix);
# my @File=glob("*edu_LGG*TCGA*494*.txt");
my @File=glob("*edu_*TCGA*.txt");
my %PancancerType=GetPanCancerType(\@File);
my %PancancerSamID=GetPanCancerSamID(\@File);
my %PancancerSamType=GetPanCancerSamType(\@File);

open ANO,$annoFile;
my @CpGList;
while(<ANO>){
my (undef,undef,undef,$cg)=split /\t/;	
push(@CpGList,$cg);
}

my (%methFileBySam,%samVersion);
foreach my $file(@File){
my @tmp=split/\.|-/,$file;
my $sampleID=join("-",@tmp[7..10]);
my $version=$tmp[4];
my $cancerType=$PancancerType{$file};
my $sampleType=substr($sampleID,13,2);
$methFileBySam{$sampleID}{$version}=$file;
push (@{$samVersion{$sampleID}},$version);
}

my %FileNewList=GetMaxVersionSamList(\@File,\%PancancerSamID);

# collect basic information from File name
open OUT1,">PancancerMethylationSampleInformation.txt";
foreach my $samid(sort keys %methFileBySam){
	my $lastVersion=$FileNewList{$samid};
  	my $fileName=$methFileBySam{$samid}{$lastVersion};
	my $cancerType=$PancancerType{$fileName};
	my $sampleType=substr($samid,13,2);
	print OUT1 "$fileName\t$samid\t$cancerType\t$sampleType\n";
}
close OUT1;


# collect methylation matrix to matrix
foreach my $file(@File){
	my @tmp=split/\.|-/,$file;
	my $sampleID=join("-",@tmp[7..10]);
	my $version=$tmp[4];
	my $cancerType=$PancancerType{$file};
	my $sampleType=substr($sampleID,13,2);
    if($sampleType eq '01' or $sampleType eq '11'){	
	open F,$file;
	while(<F>){
	chomp;
	next if /Composite/;
	next if /Hybridization/;
	my ($cg,$value,undef)=split/[\||\t]/;
	$value = sprintf("%.3f",$value) if ($value ne "NA");
	$MethDataMatrix{$file}{$cg}=$value;	
	}
	print "Start reading $file\t$cancerType\t$sampleID\t$sampleType\n";
    }
}

# write methylation matrix to all the samples by annotation files
my @FileList=sort(keys %MethDataMatrix);
my $header=join("\t",@FileList);
my $outputFile=shift @ARGV;
open OUT,">outputFile";
print OUT "\t$header\n";
foreach my $cg(@CpGList){
	print OUT "$cg";
	foreach my $file(@FileList){
	print OUT "\t$MethDataMatrix{$file}{$cg}";	
	}
	print OUT "\n";
}
close OUT;

sub GetPanCancerType{
# get cancer type by full file name		
my %PancancerCancerType;
my $fileList=shift;
my @fileList=@{$fileList};
foreach my $file(@fileList){
my @tmp=split/\.|_/,$file;
$PancancerCancerType{$file}=$tmp[2];
# print "$file\t$tmp[2]\n";
}
return(%PancancerCancerType)
}

sub GetPanCancerSamID{
# get sample id by full file name	
my %PanCancerSamID;
my $fileList=shift;
my @fileList=@{$fileList};
foreach my $file(@fileList){
my @tmp=split/\.|_|-/,$file;
$PanCancerSamID{$file}=join("-",@tmp[7..10]);
}
return(%PanCancerSamID);
}

sub GetPanCancerSamType{
# get sample type by full file name
my %PanCancerSamType;
my $fileList=shift;
my @fileList=@{$fileList};
foreach my $file(@fileList){
my @tmp=split/\.|_/,$file;
#print "$tmp[6]\t";
my $samType=substr($tmp[6],13,2);
$PanCancerSamType{$file}=$samType;
}
return(%PanCancerSamType)
}


sub GetMaxVersionSamList{
# get file list with max version No.
my (%SamMaxVersion,%MaxVersionBySamID);
my $fileList=shift;
my $PancancerSamID=shift;
my %PancancerSamID = %{$PancancerSamID};
my @fileList=@{$fileList};
foreach my $file(@fileList){
my @tmp=split/\.|-|_/,$file;
my $version=$tmp[5];
my $samID=join("-",@tmp[8..11]);
push(@{$SamMaxVersion{$samID}},$version);
# print "$samID\n";
}
foreach my $samID(sort keys %SamMaxVersion){
my @sortversion=sort {$a<=>$b} @{$SamMaxVersion{$samID}};
print "";
$MaxVersionBySamID{$samID}=$sortversion[-1];
# print "$samID\t$sortversion[-1]\n";
# print "$samID\n";
}
return(%MaxVersionBySamID);
}
