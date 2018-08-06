use strict;
use Cwd;
use POSIX;
my $dir = getcwd;
my @file=glob("*.hg38.bw");
foreach my $file(@file){
open OUT,">$file.job";
my ($fn)=split/.hg38.bw/,$file;
print OUT "#PBS -N $fn\n";
print OUT "#PBS -l nodes=1:ppn=1\n";
print OUT "cd /gpfs/home/guosa/run/blueprint\n";
print OUT "bigWigToBedGraph $file $fn.hg38.bedgraph\n";
print OUT "liftOver $fn.hg38.bedgraph /gpfs/home/guosa/hpc/db/hg38/hg38ToHg19.over.chain.gz  $fn.hg19.bedgraph\n";
print OUT "bedGraphToBigWig $fn.hg19.bedgraph /gpfs/home/guosa/hpc/db/hg19/hg19.chrom.sizes $fn.hg19.bw\n";
}
close OUT;
