
#!/usr/bin/perl
use strict;

my $fa="/home/shg047/oasis/db/hg19/hg19.fa";
my @bam=glob("CTR*bam");
my $bed=shift @ARGV;


my $output=shift @ARGV;
open F,$bed;
chomp (my @bed=<F>);
close F;

open OUT,">$output.amf";
my $header=join("\t",@bam);
print OUT "\t$header\n";
foreach my $bed(@bed){
my @tmp=split/\s+/,$bed;
my $cor="$tmp[0]:$tmp[1]-$tmp[2]";
print OUT "$cor";
foreach my $bam(@bam){
my ($sam,undef)=split /_/,$bam;  # CTR84_trimmed.fq.gz_bismark_bt2.sort.bam
my $cmd="PileOMeth extract -r $cor --minDepth 5 --mergeContext $fa $bam -o $sam.$cor";
system($cmd);
open F2,"$sam.$cor\_CpG.bedGraph" || die "Can not open the file Tmp.txt_CpG.bedGraph";
my $data;
my $count;
while(<F2>){
next if /track/i;
exit if /^\s+$/;
chomp;
my @tmp=split/\s+/;
$data+= $tmp[3];
$count++;
}
close F2;
$data =$data/($count+0.0001);
$data=sprintf("%.2f", $data);
if($count>0){
print OUT "\t$data";
}else{
print OUT "\tNA";
}
system("rm $sam.$cor\_CpG.bedGraph");
}
print OUT "\n";
}
