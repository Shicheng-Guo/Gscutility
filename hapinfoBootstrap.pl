#!/usr/bin/perl

# USAGE: perl hapinfoBootstrap.pl 100000 SRX080192.hapInfo.txt

use strict;

my $total_haps = 0;
#first count total haplotypes
open(IN, "$ARGV[1]") || die("error reading file\n");
while(my $line = <IN>){
        chomp($line);
        my @fields = split "\t", $line;
        $total_haps+=$fields[2];
}
close(IN);

my $number_of_reads_required = $ARGV[0];
#die("Too many reads requested!\n") if($number_of_reads_required > $total_haps);
my %keepPos;
my $minimum = 1;
my $maximum = $total_haps + 1;
while($number_of_reads_required > 0){
        my $x = $minimum + int(rand($maximum - $minimum));
        $keepPos{$x}++;
        $number_of_reads_required--;
}

my $cnt = 0;
open(IN, "$ARGV[1]") || die("error reading file\n");
while(my $line = <IN>){
        chomp($line);
        my $keep = 0;
        my @fields = split "\t", $line;
        while($fields[2] >0){
                $cnt++;
                if($keepPos{$cnt}){
                        $keep+=$keepPos{$cnt};
                }
                $fields[2]--;
        }
        $fields[2] = $keep;
        print join("\t", @fields), "\n" if($keep > 0);
}
close(IN);
