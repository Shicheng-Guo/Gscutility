# GSM3683988_Sample8_Swift_NovaSeq_CpG_report.txt.gz

## chr17   90      +       1       2       CG      CGA
## chr17   91      -       1       5       CG      CGT
## chr17   141     +       4       1       CG      CGA
## chr17   142     -       10      0       CG      CGC

use strict;
use Cwd;
my $file=shift @ARGV;
open F,"gunzip -c $file|" || die "cannot open the $file\n";
open OUT,">$file.bedgraph";
while(<F>){
my $c=$_;
my $g=<F>;
# print "$c";
# print "$g";
my @c=split/\s+/,$c;
my @g=split/\s+/,$g;

if($c[2]=="+" and $g[2]=="-" and $g[1]=$c[1]+1){
my $mf=($c[3]+$g[3])/($c[3]+$g[3]+$c[4]+$g[4]+0.1);
print OUT "$c[0]\t$c[1]\t$g[1]\t$mf\n";
}
}
close OUT;

