#!/bin/env perl

use strict;
use warnings;

my $usage = "Usage: $0 <bam_flag>\n";
my $flag = shift or die $usage;

die "Please enter a numerical value\n" if $flag =~ /\D+/;

if ($flag & 16){
   print "bottom strand\n";
   print "template having multiple segments in sequencing\n";
}

exit(0);

__END__
