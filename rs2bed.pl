
use strict;
use Cwd;
chdir getcwd;
my $file=shift @ARGV;
my %snp;
open F,$file || die "cannot open $file\n";
while(<F>){
my($chr,$start,$end,$rs)=split/\s+/;
my $cord="$chr:$start-$end";
$snp{$rs}=$cord;
}
close F;

my @file=glob("gwas*.txt");
foreach my $file(@file){
print "$file....\n";
open F,$file;
open OUT,">$file.bed";
while(<F>){
  chomp;
   my($snp)=split/\s+/;
   if($snp{$snp}){
   my($chr,$start,$end)=split /[:|-]/,$snp{$snp};
   print OUT "$chr\t$start\t$end\t$snp\n";
}
}
close F;
}
