#!/usr/bin/perl -w
# perl hapinfo2mf
use strict;
die USAGE() if scalar(@ARGV)<1;
my %mch_load_matrix;
my %probe_HMH_samples;
my %hap_count_matrix;
chdir shift(@ARGV);
my @hapInfo_files=glob("*.hapInfo.txt");
my @sample_list;

# read hapinfo files and save to memory 
foreach my $hapInfo_file(@hapInfo_files){
	my $sample_name = $hapInfo_file;
	$sample_name =~ s/.hapInfo.txt//;
	push(@sample_list, $sample_name);
	open(INFILE, "$hapInfo_file") || die("Error in opening $hapInfo_file!");
	while(my $line = <INFILE>){
		chop($line);
		my @fields = split(/\t/, $line);
		next if(scalar(@fields)<4);
		my $probeID = $fields[0];
		my $hapString = $fields[1];
		next if(length($hapString)<1);		
		$hap_count_matrix{$probeID}->{$sample_name}->{$hapString}=$fields[2];
	}
	close(INFILE);
}

my %aml_matrix;
my $hapStrings;
foreach my $probeID (sort keys(%hap_count_matrix)){
	foreach my $sample_name (sort keys(%{$hap_count_matrix{$probeID}})){
		my $aml=0;
                my $hapStrings;
		# calculate the average methyation level of each haplotype
		foreach my $hapString (keys(%{$hap_count_matrix{$probeID}->{$sample_name}})){
			next if($hapString =~ /[NAG]/i);
			foreach (my $i=1;$i<=$hap_count_matrix{$probeID}->{$sample_name}->{$hapString};$i++){
			$hapStrings .=$hapString;
			}
		}
		my $C_number = () = $hapStrings =~ /C/gi;
		my $T_number = () = $hapStrings =~ /T/gi;
                #print "$sample_name\t$probeID\t$C_number\t$T_number\n";
		$aml = ($C_number/($T_number+$C_number));
		$aml_matrix{$probeID}{$sample_name}=$aml;
	}
}

print "Probe_id\t", join("\t",@sample_list), "\n";
foreach my $probeID (sort keys(%aml_matrix)){
	print "$probeID";
	foreach my $sample_name(@sample_list){
		$aml_matrix{$probeID}->{$sample_name}="NA" if(! defined $aml_matrix{$probeID}->{$sample_name});
		print "\t", $aml_matrix{$probeID}->{$sample_name};
	}
	print "\n";
}

sub USAGE{
        print "\nUSAGE: perl $0 Hapinfo_Directory > output.mf\n\n";
}





