
#!/usr/bin/perl

# merge the haploinfo files created by bam2hapinfo.pl with different chrosome. 
use strict;
use Cwd;
chdir getcwd;

mkdir "../mergeHapinfo" if ! -e "../mergeHapinfo";
my @file=glob("*hapInfo.txt");

my %sam;
foreach my $file(@file){
        my ($sam,undef)=split /\./,$file;
        $sam{$sam}=$sam;
}

foreach my $sam(sort keys %sam){
my @FileTmp=glob("$sam*");
my $outFile="../mergeHapinfo/$sam.hapInfo.txt";
print "$outFile\n";
open OUT,">$outFile";
foreach my $filetmp(@FileTmp){
        open F,$filetmp;
        while(<F>){
                print OUT $_;
        }
}
close OUT;
}
