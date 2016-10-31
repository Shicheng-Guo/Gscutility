#!/usr/bin/perl -w

#Copyright 2009, 2010 The Regents of the University of California.  All Rights Reserved.

use strict;
my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
my $commonInsert = 'GTTGGAGGCTCATCGTTCCTATTCAGGCAGATGTTATCGAGGTCCGAC'; #V8 -2bp

my $AP1 = ""; # for IDT oligos
my $AP2 = ""; #

my $targetName;

while(my $line = <STDIN>){
                chop($line);chop($line);
                my ($rcH1, $rcH2)= split(/\t/, $line);
                next if(!$rcH1 || !$rcH2);
                my ($H1, $H2) = (&revComp($rcH1), &revComp($rcH2));
                my $len = length($H1) + length($H2);
                my ($chr,$chrStart,$chrEnd) = split(/[:\-]/, $exonInfo);
                my $probe = $AP1.$H1.$commonInsert.$H2.$AP2;
                print $line, "\t", $exonId, "_", $start, "\t", $probe, "\n";
}


#-------------------
sub revComp(){
#-------------------
    my $seq = shift;
        my $seqLen = length($seq);
        my $revcom = $rcTable{substr($seq,0,1)};
        for(my $i = 1; $i < $seqLen; $i++){
                $revcom = $rcTable{substr($seq, $i, 1)} . $revcom;
        }
    return $revcom;
}
