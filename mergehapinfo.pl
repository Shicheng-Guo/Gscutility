
#!/usr/bin/perl -w
use strict;
use Cwd;
use Sort::Array qw/Sort_Table/;

chdir getcwd;

my $hapfile=shift @ARGV;
open F,$hapfile;
my %data;
while(<F>){
        chomp;
        my ($mhb,$haptype,$number,$pos)=split /\t/;
        my $key="$mhb\_$haptype\_$pos";
        $data{$key}+=$number
}

open OUT,">$hapfile.SumUniq";
foreach my $key(sort keys %data){
my ($mhb,$haptype,$pos)=split /_/,$key;
print OUT "$mhb\t$haptype\t$data{$key}\t$pos\n";
}
close OUT;
