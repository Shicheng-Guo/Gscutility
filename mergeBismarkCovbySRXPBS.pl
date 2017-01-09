#!/usr/bin/perl -w

# A perl script to merge bismark cov file by SRX list.
# Contact: Shihcheng.Guo@Gmail.com
# Version 1.3
# Go to http://sra.dnanexus.com/studies/SRP028600/samples
# Select SRS and Click Related RUNS then get the Table as the input
# Run the script in the fold of coverage files created by bismark alignmetor.

use strict;
use warnings;
use Cwd;

my $file=shift @ARGV;
my %SRA;
open F,$file;
while(<F>){
chomp;
if(/(SRR\d+)/){
	my $SRR=$1;
	if(/(SRX\d+)/){
		my $SRX=$1;
		print "$SRR\t$SRX\n";
		push @{$SRA{$SRX}},$SRR;
		}
	}
}

my @cov=glob("*.cov");

my %mf;

foreach my $SRX(sort keys %SRA){
        foreach my $SRR (@{$SRA{$SRX}}){
        	foreach my $cov(@cov){
        		if($cov=~/$SRR/){
        			system("cat $cov >> $SRX.bedgraph");
        			print "$SRX combinding completed!\n";
        		}
        	}
        }
    open OUT, ">$SRX.pbs";
	print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $SRX\n";
    print OUT "#PBS -q hotel\n";
    print OUT "#PBS -l nodes=1:ppn=1\n";
    print OUT "#PBS -l walltime=168:00:00\n";
    print OUT "#PBS -o $SRX.log\n";
    print OUT "#PBS -e $SRX.err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shicheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd \$PBS_O_WORKDIR\n";
    print OUT "perl ~/bin/bedgraphUnique.pl $SRX.bedgraph\n";        
        
}

