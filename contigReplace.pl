open F,"/gpfs/home/guosa/hpc/db/hg19/hg19.contig.txt" || die print "hg19.contig.txt";
my $chr=shift @ARGV;
my %len;
while(<F>){
chomp;
if(/ID=(\d+),/){
$len{$1}=$_;
}
}
close F;

my $file="chr$chr.dose.vcf.gz";
open F1,"zcat $file |";
while(<F1>){
if(/contig=<ID=$chr/){
print "$len{$chr}\n";
}else{
print $_;
}
}

