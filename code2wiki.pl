#!/usr/bin/perl -w
# Transfer script which can be copy to wiki
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-02-21

my $file=shift @ARGV;
open F,$file;
while(<F>){
next if /^\s+$/;
chmop;
my $line=" $_";
$line=~/\s+/\s/g
print "$line\n";
}
