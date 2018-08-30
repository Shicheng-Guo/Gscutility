#!/usr/bin/perl

# transfer exm id to rs number id based on genomic coordination.

use strict;
use Cwd;
my $dir = getcwd;
chdir $dir;

my $bim="/home/local/MFLDCLIN/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.bim";
my $allsnp=" /home/local/MFLDCLIN/guosa/hpc/db/hg19/allsnp150.hg19";

open F1,"$allsnp";
open F2,"$bim";
open OUT,">exm2rs.rlt.txt";
my %bim;
my %allsnp;

#open and read all snps from 1000 genome database.
while(<F1>){
        chomp;
        my($chr,undef,$pos,$rs,undef)=split/\s+/;
    $allsnp{"$chr:$pos"}=$rs;
}
print "SNP file reading completed.....\n";

#open bim files and read probes of immumina exom array
while(<F2>){
        chomp;
        my $line=$_;
        my($chr,$exmid,undef,$end,$ref,$alt)=split/\s+/,$line;
        my $loc="chr$chr:$end";
        if ($allsnp{$end}){
        print OUT "$chr\t$allsnp{$loc}\t0\t$end\t$ref\t$alt\n";
        }else{
        print OUT "$line\n";
        }
}
print "change the exm id to rs number completed...\n";
