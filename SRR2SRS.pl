
#!/usr/bin/perl -w

# A perl script to build config file for SRS Bam Merge
# Contact: Shihcheng.Guo@Gmail.com
# Version 1.3
# Go to http://sra.dnanexus.com/studies/SRP028600/samples
# Select SRS and Click Related RUNS then get the Table as the input

use strict;
use warnings;
use Cwd;

my $file=shift @ARGV;
my %SRA;
open F,$file;
while(<F>){
chomp;
my ($SRR,$SRP,$SRX,$SRS)=split /\t/;
push @{$SRA{$SRS}},$SRR;
}

foreach my $SRS(sort keys %SRA){
        print "$SRS";
        foreach my $SRR (@{$SRA{$SRS}}){
                print "\t$SRR";
        }
        print "\n";
}
