#!/usr/bin/perl -w

# check which files are not succeed in bam2methyfreq
# Run the script to the bam directory (err file of PBS)
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-03-23

use strict;
use Cwd;
my $dir=getcwd;

my @err=glob("*bam2mf.err");

foreach my $err(@err){
my $suc='Fail';
open F,$err;
while(<F>){
if (/Finished BedGraph conversion/){
$suc='Succeed';
last;
}
}
if($suc eq 'Fail'){
print "$err\n";
}
}
