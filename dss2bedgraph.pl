use strict;
# read all the dss files of same group and return the average methylation level
use Cwd;
my $dir=getcwd;
my @file=glob("*OA*.dss");
my %data;
foreach my $file (@file){
open F,$file;
while(<F>){
        chomp;
    my($chr,$pos,$X,$M)=split/\s+/;
    $data{"$chr:$pos"}{"X"} +=$X;
    $data{"$chr:$pos"}{"M"} +=$M;
}
}

foreach my $pos(sort keys %data){
        my($chr,$end)=split /\:/,$pos;
        my $start=$end-1;
        print "$chr\t$start\t$end\t".$data{$pos}{"M"}/$data{$pos}{"X"}."\n" if $data{$pos}{"X"}>4;
}
