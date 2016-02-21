#!/usr/bin/perl -w

# estimate the length of the Introns from GENCODE GTF file
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-02-21
# Contact: Shicheng Guo<scguo@ucsd.edu>; Kun Zhang<kzhang@eng.ucsd.edu>
use strict;
use Cwd;
chdir getcwd;

my $input=@ARGV[0];

open F,$input;
print "\{|{{table}}\n";
while(<F>){
chomp;
my @line=split /\t/;
if($. eq 1){
print "{|style=\"font-size:80%;\"\m"
foreach my $ele(@line){
print "|align=\"center\" style=\"background:#f0f0f0;\"|\'\'\' $ele\'\'\'\n";	
}
print "|-\n";
}else{
my $tmp=join("||",@line);
print "| $tmp\n|-\n";
}
print"|\}\n";
}
