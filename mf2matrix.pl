#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;
# usage: bedGraph2matrix.pl >  bedGraph_matrix.txt
# collect all the bedGraph to matrix
print "\n Usage: perl bedGraph2matrix.pl ./\n";
print "The output file will be printed to project.rlt.txt\n";


my @bedgraph_files=glob("*.bw.mf");
my %sample;
foreach my $file(@bedgraph_files){
        my ($sample,undef)=split /\./,$file;
        $sample{$sample}=$sample;
}

my %mf_matrix;
my %coor;
foreach my $sample(sort keys %sample){
        my @file=glob("$sample*.bw.mf");
        foreach my $file(@file){
                open F,$file;
                while(<F>){
			chomp;
                        next if /track/;
                        my ($coor,$len,$cov,$avg1,$avg2,$mf)=split /\s+/;
                        if($cov>0){
                        $mf_matrix{$sample}{$coor}=$mf;
                        }else{
                      	$mf_matrix{$sample}{$coor}="NA";	
                        }
                        $coor{$coor}=$coor;
                }
        }
}

my $head=join("\t",sort keys %sample);
print "\t$head\n";
        foreach my $coor(sort keys %coor){
                print "$coor";
                foreach my $sample(sort keys %sample){
                        if(defined $mf_matrix{$sample}{$coor}){
                        print "\t$mf_matrix{$sample}{$coor}";
                        }else{
                        print "\tNA";
                }
        }
        print "\n";
}

