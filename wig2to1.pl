use strict;
my $file=shift @ARGV;
open F,$file || die "cannot find $file\n";
my %cm;
my %ct;
my %pos;
my $chrome;
while(<F>){
  next if /CHROMOSOME/;
  chomp (my $line=$_);
  my ($chr,$pos,$ref,undef,undef,undef,undef,$cm,undef,$cu)=split/\s+/,$line;
  if($ref eq "C"){
  $cm{$chr}{$pos}+=$cm;
  $ct{$chr}{$pos}+=($cm+$cu);
  $pos{$chr}{$pos}=$pos;
  }else{
  $cm{$chr}{$pos-1}+=$cm;
  $ct{$chr}{$pos-1}+=($cm+$cu);
  $pos{$chr}{$pos-1}=$pos-1;
  }
}

foreach my $chr(sort keys %pos){
        foreach my $pos(sort keys %{$pos{$chr}}){
        print "$chr\t$pos\t$ct{$chr}{$pos}\t$cm{$chr}{$pos}\n";
        }
}
