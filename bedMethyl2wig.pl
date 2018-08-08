use strict;
my $input=shift @ARGV;
open F,$input;
open OUT,">$input.bedgraph";
while(<F>){
        chomp(my $F1=$_);
        chomp(my $F2=<F>);
        my($chr,$start,undef,undef,undef,undef,undef,undef,undef,$read1,$mf1)=split/\s+/,$F1;
        my(undef,$check,$end,undef,undef,undef,undef,undef,undef,$read2,$mf2)=split/\s+/,$F2;
        my $readsum=$read1+$read2;
        print "$start\t$check\tError,CpG dinucleotide assumption was broken!\n" if $check !=$start+1;
        next if $readsum<3;
        my $mf=($read1*$mf1/100+$read2*$mf2/100)/$readsum;
        print OUT "$chr\t$start\t$end\t$mf\n";
}
