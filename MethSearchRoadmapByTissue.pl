#!/usr/bin/perl -w
use strict;
use Cwd;
use Statistics::Basic qw(:all);
# Search Methylation Status from Roadmap
# Contact: Shicheng Guo(Shihcheng.Guo@gmail.com)
# Esophagus=23
# Version 1.3
# Update: 2016-05-12

die &USAGE if scalar @ARGV<4;

my $column=23;                   # Esophagus
my $gap=500;

my $input=shift @ARGV;
my $OUTPUT=shift @ARGV;
$column=shift @ARGV;
$gap=shift @ARGV;

open F,$input;
my @cor;
while(<F>){
my ($chr,$start,$end)=split/\s+/;	
my $cor="$chr:$start-$end";
push(@cor,$cor);
}

open OUT,">$OUTPUT.fm.rlt.txt";
foreach my $cor(@cor){
my ($chr,$start,$end)=split /:|-/,$cor;
open F2,"$chr.fm" || die "cannot open $chr.fm";
my $num;
my $sum;
my @fm;
while(<F2>){
chomp;
my @line=split/\s+/;
if($line[0]>=$start-$gap && $line[0]<=$end+$gap){
$num++;
push(@fm,$line[$column]);
print OUT "$chr\t$start\t$end\t$line[0]\t$line[$column]\n";	
}
if($line[0]>$end){
last;
}	
}
my $avg = median(@fm);
my $sd = stddev(@fm);
my $cpgnum=int($num/2);
print "$chr\t$start\t$end\t$cpgnum\t$avg\t$sd\n";
close F2;
}

sub USAGE{
print "\nUsage: perl $0 TargetBed Prefix ColumnID GAP\n";
print "Search Methylation Status from Roadmap\n";

print '
Format For: TargetBed: 
chr12   95942750        95942970
chr17   43339200        43339400
chr6    110678910       110679110
chr6    133562350       133562550
chr8    67344530        67344730
chr16   51184245        51184465
';

print '
Roadmap Dataset Download: 
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/EG.mnemonics.name.xls
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/FractionalMethylation.tar.gz
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/FractionalMethylation.tar.gz.md5sum
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/header
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/ReadCoverage.tar.gz
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/ReadCoverage.tar.gz.md5sum
';

print '
FM file Header: 
E003    ESC.H1  H1_Cell_Line
E004    ESDR.H1.BMP4.MESO       H1_BMP4_Derived_Mesendoderm_Cultured_Cells
E005    ESDR.H1.BMP4.TROP       H1_BMP4_Derived_Trophoblast_Cultured_Cells
E006    ESDR.H1.MSC     H1_Derived_Mesenchymal_Stem_Cells
E007    ESDR.H1.NEUR.PROG       H1_Derived_Neuronal_Progenitor_Cultured_Cells
E008    ESC.H9  H9_Cell_Line
E011    ESDR.CD184.ENDO hESC_Derived_CD184+_Endoderm_Cultured_Cells
E012    ESDR.CD56.ECTO  hESC_Derived_CD56+_Ectoderm_Cultured_Cells
E013    ESDR.CD56.MESO  hESC_Derived_CD56+_Mesoderm_Cultured_Cells
E016    ESC.HUES64      HUES64_Cell_Line
E017    LNG.IMR90       IMR90_Cell_Line
E021    IPSC.DF.6.9     iPS_DF_6.9_Cell_Line
E022    IPSC.DF.19.11   iPS_DF_19.11_Cell_Line
E024    ESC.4STAR       4star
E050    BLD.MOB.CD34.PC.F       Mobilized_CD34_Primary_Cells_Female
E053    BRN.CRTX.DR.NRSPHR      Neurosphere_Cultured_Cells_Cortex_Derived
E054    BRN.GANGEM.DR.NRSPHR    Neurosphere_Cultured_Cells_Ganglionic_Eminence_Derived
E058    SKIN.PEN.FRSK.KER.03    Penis_Foreskin_Keratinocyte_Primary_Cells_skin03
E065    VAS.AOR Aorta
E066    LIV.ADLT        Adult_Liver
E070    BRN.GRM.MTRX    Brain_Germinal_Matrix
E071    BRN.HIPP.MID    Brain_Hippocampus_Middle
E079    GI.ESO  Esophagus
E084    GI.L.INT.FET    Fetal_Intestine_Large
E085    GI.S.INT.FET    Fetal_Intestine_Small
E094    GI.STMC.GAST    Gastric
E095    HRT.VENT.L      Left_Ventricle
E096    LNG     Lung
E097    OVRY    Ovary
E098    PANC    Pancreas
E100    MUS.PSOAS       Psoas_Muscle
E104    HRT.ATR.R       Right_Atrium
E105    HRT.VNT.R       Right_Ventricle
E106    GI.CLN.SIG      Sigmoid_Colon
E109    GI.S.INT        Small_Intestine
E112    THYM    Thymus
E113    SPLN    Spleen
';

}


