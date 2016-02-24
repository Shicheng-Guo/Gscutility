#!/usr/bin/perl
#table2wikitable.pl
use strict;
use Cwd;
chdir getcwd;

my $input=@ARGV[0];

open F,$input;
while(<F>){
chomp;
my @line=split /\t/;
if($. eq 1){
print "{|style=\"font-size:80%;\"\n";
foreach my $ele(@line){
print "|align=\"center\" style=\"background:#f0f0f0;\"|\'\'\' $ele\'\'\'\n";
}
print "|-\n";
}else{
my $tmp=join("||",@line);
print "| $tmp\n|-\n";
}
}

print"|\}\n";
