#!/usr/bin/env perl
use strict;
use warnings;
use Bib::Tools qw(:all);
my $refs = Bib::Tools->new();
$refs->add_google_search('https://scholar.google.com/scholar?as_ylo=2019&q=rheumatoid+arthritis&hl=en&as_sdt=0,50');
print $refs->print;
