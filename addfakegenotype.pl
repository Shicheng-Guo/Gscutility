use strict;
use Cwd;
chdir getcwd;
open F,"dbSNP153.hg19.vcf";
while(<F>){
chomp;
my $line=$_;
if(/^#CHROM/){
print "$_\tFORMAT\tGSC001\n" if /^#/;
}elsif(/^#/){
print "$_\n" if /^#/;
}else{
$line=~s/NC_(0)+//i;
$line=~s/\.\d+//i;
print "$line\tGT\t0|1\n";
}
}
