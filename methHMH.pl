#!/usr/bin/perl -w
use strict;
use Cwd;
chdir getcwd;

# Extract Pileup within specfic genomic regions (CDS, Exon)
# Set Pileup and GenomicInterval regions
# Contact: Shicheng Guo
# Version 1.3
# Update: 2017-04-5

my %sam;
my %hapinfo;
my %loc;
my %hap;
my %cpg;

 die &USAGE if scalar @ARGV<3;
 my $input=shift @ARGV;
 my $bed=shift @ARGV;
 my $output=shift @ARGV;
# my $input="input.txt";
# my $output="output.txt";

open F,$bed;
my $window=50;
my %hash;

while(<F>){
        chomp;
        next if /^\s+$/;
        my ($chr,$start,$end,$id,$gene,$block)=split/\s+/;
        my $bin1=int($start/$window);
        my $bin2=int($end/$window);
        foreach my $i($bin1..$bin2){
        $hash{$chr}{$i-1}="$chr:$start-$end";
        $hash{$chr}{$i}="$chr:$start-$end";
        $hash{$chr}{$i+1}="$chr:$start-$end";
		# print "$i\n";
	}
}
close F;

open F1,$input || die "cannot open $input\n";
print "\nStart to reading haplotype files....\n";
while(<F1>){
    chomp;
    next if /^\s*$/;
    my ($sam,$file)=split/\t/;
    print "\t$file reading completed....\n";
	open F2,$file || die "cannot open $file\n";
	while(<F2>){
	chomp;
	my ($loc,$haptype,$count,$cpg)=split /\s+/;
    next if length($haptype)<2;    # TBD
	my ($chr,undef,undef)=split/[:-]/,$loc;
	my ($start,undef)=split/,/,$cpg;
	my $bin=int($start/$window);
    if(defined $hash{$chr}{$bin}){
	    $hapinfo{$loc}{$cpg}{$haptype}{$sam}+=$count;
		$loc{$loc}=$loc;
		$cpg{$cpg}=$cpg;
		$hap{$haptype}=$haptype;
		$sam{$sam}=$sam;
	}
	}
	close F2;
}
close F1;

# print  
open OUT1,">$output.cctf";
open OUT2,">$output.homf";
my @sam=sort keys %sam;
my $header=join("\t",@sam);
# print OUT1 "GenomeInterval\tCpGLocList\t$header\n";
# print OUT2 "GenomeInterval\tCpGLocList\t$header\n";
print OUT1 "Chr\tStart\tEnd\t$header\n";
print OUT2 "Chr\tStart\tEnd\t$header\n";

print "Start to Performance Haplotype Split and Recombination....\n";

my @loc=keys %loc;
my $loci=0;
my $part=(int($#loc/10+1));
# print "$#loc\tpart=$part\n";
foreach my $loc(sort keys %hapinfo){
	## Step 1. performance haplotype split and recombination
    $loci++;	
	my %Hapinfo;
    foreach my $cpg(sort keys %{$hapinfo{$loc}}){
	foreach my $hap(sort keys %{$hapinfo{$loc}{$cpg}}){
		my @cpg=split/,/,$cpg;
		my $i=0;
		while($i<($#cpg)){
			my $j=2;
			while($j<=$#cpg-$i+1){
			foreach my $sam(sort keys %sam){
			my $string_hap=substr($hap,$i,$j);
			my $string_cpg=join(",",@cpg[$i..($i+$j-1)]);
			$hapinfo{$loc}{$string_cpg}{$string_hap}{$sam}=$hapinfo{$loc}{$cpg}{$hap}{$sam};
			if(defined $hapinfo{$loc}{$cpg}{$hap}{$sam}){
			$Hapinfo{$loc}{$string_cpg}{$string_hap}{$sam}+=$hapinfo{$loc}{$string_cpg}{$string_hap}{$sam};
			}
			}
			$j++;
			}
			$i++;
			}
		#print "\n";
		}
	}
	
	## Step 2. summary and report homozygote haplotype freqeunce.
	foreach my $loc(sort keys %Hapinfo){
        # print "$loc\n";
        my ($chr,undef)=split/:/,$loc;
		foreach my $cpg(sort keys %{$Hapinfo{$loc}}){
            my @tmp=split/,/,$cpg;
            my $bin=int($tmp[0]/$window);
			next if ! defined $hash{$chr}{$bin};
			my (undef,$start,$end)=split/[:|-]/,$hash{$chr}{$bin};
			print OUT1 "$chr\t$start\t$end";
			print OUT2 "$chr\t$start\t$end";
			foreach my $sam(sort keys %sam){
				my $C=0;
				my $T=0;
				my $N=0;
				my $ratio1;
				my $ratio2;
				foreach my $hap(sort keys %{$Hapinfo{$loc}{$cpg}}){
                next if length($hap)<2;
				if (! defined $Hapinfo{$loc}{$cpg}{$hap}{$sam}){
				$ratio1="NA";
				$ratio2="NA";
				}else{	
				my $num=length($hap);
    			if($hap eq "C"x $num){
    				$C=$C+$Hapinfo{$loc}{$cpg}{$hap}{$sam};
    			}elsif($hap eq "T"x $num){
    				$T=$T+$Hapinfo{$loc}{$cpg}{$hap}{$sam};
    			}else{
					$N=$N+$Hapinfo{$loc}{$cpg}{$hap}{$sam};
				}
				$ratio1=sprintf("%.2f",$C/($T+$C+0.001));
				$ratio2=sprintf("%.2f",(1-($N/($T+$C+$N+0.001))));
				}
				}
				print OUT1 "\t$ratio1";
				print OUT2 "\t$ratio2";
			}
			print OUT1 "\n";
			print OUT2 "\n";

		}
	}
	
	# print "$loci\n";
	my $status=sprintf("\t........%2d%%(%d/%d) completed!",100*$loci/($#loc+1),$loci,$#loc+1);
	print "$status\n" if $loci % $part eq 0;
}
print "\n";	
close OUT1;
close OUT2;



sub USAGE{
print "Usage: perl $0 samplesheet.txt interest.bed prefix\n";
print "Methylation Haplotype and Counts in Each File will be Collected with Mode of Matrix\n";

print '
Format For: samplesheet file: 
Indx01  Indx01.sortc.bam.hapInfo.txt
Indx02  Indx02.sortc.bam.hapInfo.txt
';
print '
Format For: Interest genomic interval file: 
chr10   76532564        76532591        Col18a1 Exon
chr10   76532243        76532564        Col18a1 Intron
';
}




