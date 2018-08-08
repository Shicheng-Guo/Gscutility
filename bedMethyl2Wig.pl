use strict;
my $input=shift @ARGV;
open F,$input;
open OUT,">$input.bedgraph";
while(<F>){
        chomp;
        my($chr,$start,$end,undef,undef,undef,undef,undef,undef,$read,$mf)=split/\s+/;
        next if $read<1;
        $mf=$mf/100;
        print OUT "$chr\t$start\t$end\t$mf\n";
}
