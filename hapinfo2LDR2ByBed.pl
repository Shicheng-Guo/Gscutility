#!/usr/bin/perl
use Cwd;
chdir getcwd;
printUsage() if scalar(@ARGV)<2;
my $bed=shift @ARGV;
my $hapinfo=shift @ARGV;
open F1,$bed;
while(<F1>){
chomp;
next if /^\s+$/;
my ($chr,$start,$end)=split/\t/;
my $cor="$chr:$start-$end";
system("perl hapinfo2LDR2.pl tmp $cor \< $hapinfo");
open F2,"tmp.$chr.rsq";
my $i=0;
my $Rlt=0;
# print "----------------------------------------------------------\n";
while(<F2>){
# print "$_";
next if /^\t/;
my (undef,@line)=split/\s+/;
foreach my $tmp(@line){
next if $tmp=~/NA/i;
#print("$i\t$tmp\t$Rlt\n");
$i++;
$Rlt=$Rlt+$tmp;
}
}
# print "-----------------------------------------------------------\n";
my $output=$Rlt/$i;
print "$cor\t$output\n";
}


sub printUsage{
        print " \nUsage:\n";
                print " perl $0 input.bed input.hapinfo\n";
                print " For example:\n perl $0 chr1.bed haploinfo.input.txt\n\n";
                print "-----------------------------------------------------------------------------\n";
                print " haploinfo.input.txt format:\n";
                print " chr10:10000873-10001472        CCC     1       10001056,10001082,10001168\n\n";
        print "-----------------------------------------------------------------------------\n";
        print " Output Format:\n";
        print " chr10:10000873-10001472\t0.64\n\n";
        print "-----------------------------------------------------------------------------\n";

        exit 0;
}
