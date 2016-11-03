#!/usr/bin/perl -w
use strict;
use Cwd;

# Build layered UCSC tracks to bigWig files in 
# Set Pileup and GenomicInterval regions
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-03-16
my @file=glob("*bw");
foreach my $file(@file){
my ($id,$name)=split /\./,$file;
my $iid=$1 if ($id=~/(GSM\d+)/);

next if $file=~/Mobilized/i;  # Mobilized
next if $file=~/GSM1279520/i; # Metastasis
next if $file=~/U87MG/i; # U87MG is mis-lable in ATCC
next if $file=~/H1_Derived/i; # U87MG is mis-lable in ATCC
next if $file=~/H1-BMP4/i; # ESC differentiation

# General
$name=~s/_\d+//ig;
$name=~s/Left_//ig;
$name=~s/Right_//ig;
$name=~s/_Tissue//ig;

# Stomach
$name=~s/Gastric_/Stomach/ig;

# Fetal or Adult
$name=~s/Fetal_//ig;
$name=~s/Adult_//ig;

# Colon
$name=~s/Colon_P/Colon-Cancer/ig;
$name=~s/Sigmoid_Colon/Colon/ig;

# Brain
$name=~s/Brain_Hippocampus_Middle/Brain-Normal/ig;
$name=~s/Brain_Germinal_Matrix/Brain-Normal/ig;
$name=~s/Brain_W/Brain-Normal/ig;

# Breast
$name=~s/Breast_Luminal_Epithelial_Cells/Breast/ig;
$name=~s/Breast_Myoepithelial_Cells/Breast/ig;
$name=~s/Breast/Breast-Normal/ig;

# Lung
$name=~s/H1437/LUAD/ig;
$name=~s/H1672/SCLC/ig;

# Oral 
$name=~s/H157/OSCC/ig;

# ESC
$name=~s/HUES64/ESC/ig;
$name=~s/iPS_DF/iPS/ig;

# Muscle
$name=~s/Psoas_Muscle/Muscle/ig;
$name=~s/Muscle_Leg/Muscle/ig;

#Testis
$name=~s/Testis_Spermatozoa_Primary_Cells/Testis/ig;

# Foreskin
$name=~s/Penis_Foreskin_Fibroblast_Primary_Cells/Foreskin/ig;
$name=~s/Penis_Foreskin_Keratinocyte_Primary_Cells/Foreskin/ig;


my $newName="$iid.$name.bw";
print "$file\t$iid\t$name\t$newName\n";

}
