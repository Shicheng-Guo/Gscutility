#!/usr/bin/perl -w

# merge the list of bam files together, using one of the files's header
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-01-30
use strict;

my $out_bam = $ARGV[0];
my $cmd = "samtools merge -h";
my $line = <STDIN>;
chomp($line);
$cmd = $cmd . " $line $out_bam";
while($line = <STDIN>){
        chomp($line);
        $cmd = $cmd . " $line";
}
print $cmd, "\n";
