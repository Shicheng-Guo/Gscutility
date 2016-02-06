#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;

my @file=glob("*fq*");
my %data;
foreach my $file(@file){
my ($sample,undef)=split /\./,$file;
my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks)= stat($file);
if ($file=~/gz/){
$data{$sample}{'gz'}=$size;
}else{
$data{$sample}{'ngz'}=$size;
}
}

open OUT,">compress.ratio.txt";
foreach my $sample(sort keys %data){
        print OUT "$sample\t";
        my $ratio=$data{$sample}{'gz'}/$data{$sample}{'ngz'};
        print OUT "$data{$sample}{'gz'}\t$data{$sample}{'ngz'}\t$ratio";
        print OUT "\n";
}
