# How to add fake genotypes for dbsnp152 vcf file

use strict;
open F,"/home/guosa/hpc/db/hg19/dbSNP152.hg19.vcf";
while(<F>){
chomp;
my $line=$_;
if(/^#CHROM/){
print "$_\tGSC001\n" if /^#/;
}elsif(/^#/){
print "$_\n" if /^#/;
}else{
$line=~s/NC_000001.//i;
print "$line\t0|1\n";
}
}
