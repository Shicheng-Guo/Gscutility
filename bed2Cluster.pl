#!/usr/bin/perl -w
# collect cluster from bed files.
# Contact: Dr. Kun Zhang
# Version 1.3
# Update: 2016-02-18
# Written by Kun Zhang (kzhang@bioeng.ucsd.edu), last modified 10/07/2016 

my $BED_file = $ARGV[0];
my $min_sites_per_cluster=4;
my $cluster_max_CpG_spacing = 100;
open(INFILE, "$BED_file")||next;
my ($curr_chr, $curr_start, $last_pos, $cluster_id,$N_sites)=(0,0,0,1,0);
while(my $line = <INFILE>){
	chop($line);
	my @fields = split(/\t/, $line);
	my $chr=$fields[0];
	my $pos=$fields[1];
	if(!$curr_chr || !$curr_start || $curr_chr ne $chr){
		$curr_chr = $chr;
		$curr_start=$pos;
		$last_pos = $pos;
		$N_sites=1;
		next;
	}elsif($curr_chr eq $chr && $pos-$last_pos < $cluster_max_CpG_spacing){
		$N_sites++;
		$last_pos = $pos;
	}else{
		print "$chr\t$curr_start\t$last_pos\t$N_sites\t$cluster_id\n" if($N_sites>$min_sites_per_cluster && $last_pos-$curr_start>60 && $chr =~ /chr[1-9]+/);
		$N_sites=1;
		$last_pos = $pos;
		$curr_start=$pos;
		$cluster_id++;
	}
}
