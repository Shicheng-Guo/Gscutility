#!/usr/bin/perl -w

# Version update: April 5th 2016 

use strict;

my $samtools = "samtools";
my $target_list_file= $ARGV[0];
my $merged_bam_file = $ARGV[1];
my $Aligner=$ARGV[2];
my $chrsizeFile=$ARGV[3];
my $cpg_position_file=$ARGV[4];
die("USAGE: mergedBam2hapInfo.pl InputBed MergeBam <bisReadMapper|bismark> ChrosmeZizeFile CpG_Position_File\n") if(scalar(@ARGV)<5);
my $phred=33;

my $target_flanking_len = 80;
my %targetTable;
my %hapInfoTable;
my %binnedtargetTable;
my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'R'}='Y';
$rcTable{'Y'}='R';
$rcTable{'M'}='K';
$rcTable{'K'}='M';
$rcTable{'S'}='S';
$rcTable{'W'}='W';
my ($reads_parsed, $reads_after_filtering);

open F,"$chrsizeFile" || die "can not open $chrsizeFile,please make sure this file is existed!";
my %chrSizes;
while(<F>){
chomp;
my($chr,$len)=split/\s+/;
$chrSizes{$chr}=$len;
}

open(INFILE, "$target_list_file")||die("Error in opening $target_list_file\n");
#forward strand
my $bin_size = 10000;
while(my $line = <INFILE>){
	chop($line);
        next if $line=~/^\s+$/;
	my ($chr,$target_start,$target_end,$id) = split(/[\t ]+/, $line);
	my $index = int($target_start/$bin_size);
	my $target_id = "$chr:$target_start-$target_end";
	push(@{$binnedtargetTable{$chr}->{$index}},$target_id);
	push(@{$binnedtargetTable{$chr}->{$index-1}},$target_id);
	push(@{$binnedtargetTable{$chr}->{$index+1}},$target_id);
}
close(INFILE);

open(INFILE, "$cpg_position_file")|| die "Error in opening $cpg_position_file, Please make sure this file existed!\n";

while(my $line = <INFILE>){
	chop($line);
	my ($chr, $pos) = split(/[\t ]+/, $line);
	my $index = int($pos/$bin_size);
	next if(!$binnedtargetTable{$chr}->{$index});	
	foreach my $target_id (@{$binnedtargetTable{$chr}->{$index}}){
		my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
		if($target_start <= $pos && $pos <=$target_end){
			push(@{$targetTable{$target_id}->{"CpG_positions"}},$pos);			
		}
	}
}
close(INFILE);

if($Aligner eq 'bismark'){
foreach my $target_id (sort keys(%targetTable)){
	my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
	my $region_start = $target_start-$target_flanking_len;
	my $region_end = $target_end+$target_flanking_len;
	my $cmd = "$samtools view -q 10 $merged_bam_file $chr:$region_start-$region_end";
	my $ans = `$cmd`;
	my @lines = split(/\n/, $ans);
	next if(scalar(@lines)<1 || !$targetTable{$target_id}->{"CpG_positions"});
	$reads_parsed += scalar(@lines);
	$reads_after_filtering += scalar(@lines);
	my @CpG_positions = sort {$a <=> $b} @{$targetTable{$target_id}->{"CpG_positions"}};
	my %readHapInfo;	

	#print "##$cmd\n";
	#Going through the reads one at a time
	my %read_start_pos_table;
	foreach my $line (@lines){
		my @fields = split(/\t/, $line);
		my $read_strand = $fields[1];		
		my $seq = $fields[9];		
		my $qual_string = $fields[10];
		my $read_start = $fields[3];
		my $CIGAR = $fields[5];
		my $read_length = length($seq);
		next if($CIGAR =~ /[ID]/);
		if($CIGAR=~ /S/){
			my ($clip_len, @others) = split(/S/, $CIGAR);
			next if(length($clip_len)>2);
			$seq=substr($seq,$clip_len,$read_length-$clip_len);		
			$qual_string=substr($qual_string,$clip_len,$read_length-$clip_len);	
			$read_length -= $clip_len;			
		}elsif($seq =~ /^CG/){
			$seq=substr($seq,2,$read_length-2);		
			$qual_string=substr($qual_string,2,$read_length-2);	
			$read_start += 2;
			$read_length -= 2;
		}
		my @read_fields = split(/[\:\#\_]+/, $fields[0]);
		my $N_read_fields = scalar(@read_fields);
		for(my $i=$N_read_fields; $i<=4; $i++){
			push(@read_fields,"NA");
		}
		my $read_id = scalar(@read_fields)>7 ? join("-", @read_fields[0..6]) : join("-", @read_fields[0..4]);	
		#print "#$read_id\t$CIGAR\t$read_start\t$read_strand\t$seq\t$qual_string\n";		
		foreach my $CpG_position(@CpG_positions){
			my $offset = $CpG_position-$read_start+1;
			next if($offset<0 || $offset >= length($seq));	
            # the situation sometimes should be change dependent on Flag of different alignmentor
			if($read_strand <16 || $read_strand ==99 ||$read_strand ==147)
            {
                                # positive chain
				my $qual_score = ord(substr($qual_string,$offset,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score); # UMIs choose the best one
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=substr($seq,$offset,1);
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;
			}else{
				# Negative chain
				my $qual_score = ord(substr($qual_string,$offset+1,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score); # UMIs choose the best one
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=$rcTable{substr($seq,$offset+1,1)};			
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;	

			}	
			#print "$offset\t", $readHapInfo{$read_id}->{$CpG_position}, "\n";
		}	
		$read_start_pos_table{$read_id} = $read_start if(!$read_start_pos_table{$read_id} || $read_start_pos_table{$read_id}<$read_start);
	}		
	
	#Use read start as the UMI to collapse multiple clonal reads.
	my @read_list = keys(%readHapInfo);		
	my %unique_read_base_info;
	foreach my $read_id (@read_list){		
		my $UMI = $read_start_pos_table{$read_id};
		foreach my $CpG_position (keys(%{$readHapInfo{$read_id}})){
			next if(!$readHapInfo{$read_id}->{$CpG_position}->{"base"});
			#print "$CpG_position\n";
			#print $readHapInfo{$read_id}->{$CpG_position}->{"base"}, "\n";
			#print $readHapInfo{$read_id}->{$CpG_position}->{"qual"}, "\n";
			$unique_read_base_info{$UMI}->{$CpG_position}->{$readHapInfo{$read_id}->{$CpG_position}->{"base"}}+=$readHapInfo{$read_id}->{$CpG_position}->{"qual"};
		}
	}
	
	#Derive the consensus haplotype string based on multiple clonal reads.
	foreach my $UMI (keys(%unique_read_base_info)){
		my $hap_string="";
		foreach my $CpG_position(@CpG_positions){
			my $base = "N";
			my $best_qual = 0;
			if($unique_read_base_info{$UMI}->{$CpG_position}){
				foreach my $base_call (keys(%{$unique_read_base_info{$UMI}->{$CpG_position}})){
					next if($unique_read_base_info{$UMI}->{$CpG_position}->{$base_call} <= $best_qual);
					$best_qual = $unique_read_base_info{$UMI}->{$CpG_position}->{$base_call};
					$base = $base_call;
				}
			}					
			$hap_string=$hap_string.$base;
		}	
		$hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string}++ if($hap_string);
	}

	#report haplotype strings only on the positions with valid calls.
	foreach my $hap_string(keys %{$hapInfoTable{$target_id}->{"hap_counts"}}){
		my @valid_positions;
		my $valid_hap="";
		for(my $i=0; $i<length($hap_string); $i++){
			my $allele = substr($hap_string,$i,1);
			next if($allele eq "N");
			$valid_hap = $valid_hap.$allele;
			push(@valid_positions,$CpG_positions[$i]);
		}
		next if(scalar(@valid_positions)<1);
		next if($valid_hap =~ /[AG]/);
		print $target_id,"\t$valid_hap\t", $hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string},"\t",join(",",@valid_positions), "\n";
	}
}
}else{
foreach my $target_id (sort keys(%targetTable)){
	my ($chr,$target_start,$target_end) = split(/[:\-]/, $target_id);
	my $region_start = $target_start-$target_flanking_len;
	my $region_end = $target_end+$target_flanking_len;
	my $cmd = "$samtools view -q 10 $merged_bam_file $chr:$region_start-$region_end";
	my $ans = `$cmd`;
	my @lines = split(/\n/, $ans);
	next if(scalar(@lines)<1 || !$targetTable{$target_id}->{"CpG_positions"});
	$reads_parsed += scalar(@lines);
	$reads_after_filtering += scalar(@lines);
	my @CpG_positions = sort {$a <=> $b} @{$targetTable{$target_id}->{"CpG_positions"}};
	my %readHapInfo;	

	#print "##$cmd\n";
	#Going through the reads one at a time
	my %read_start_pos_table;
	foreach my $line (@lines){
                #D print "$line\n";
		my @fields = split(/\t/, $line);
		my $read_strand = $fields[1];		
		my $seq = $fields[9];		
		my $qual_string = $fields[10];
		my $read_start = $fields[3];
		my $CIGAR = $fields[5];
		my $read_length = length($seq);
		next if($CIGAR =~ /[ID]/);
		if($CIGAR=~ /S/){
			my ($clip_len, @others) = split(/S/, $CIGAR);
			next if(length($clip_len)>2);
			$seq=substr($seq,$clip_len,$read_length-$clip_len);		
			$qual_string=substr($qual_string,$clip_len,$read_length-$clip_len);	
			$read_length -= $clip_len;			
		}elsif($seq =~ /^CG/){
			$seq=substr($seq,2,$read_length-2);		
			$qual_string=substr($qual_string,2,$read_length-2);	
			$read_start += 2;
			$read_length -= 2;
		}
		my @read_fields = split(/[\:\#\_]+/, $fields[0]);
		my $N_read_fields = scalar(@read_fields);
		for(my $i=$N_read_fields; $i<=4; $i++){
			push(@read_fields,"NA");
		}
		my $read_id = scalar(@read_fields)>7 ? join("-", @read_fields[0..6]) : join("-", @read_fields[0..4]);	
		# print "#$read_id\t$CIGAR\t$read_start\t$read_strand\t$seq\t$qual_string\n";		
		foreach my $CpG_position(@CpG_positions){
			my $offset = $CpG_position-$read_start+1;
			next if($offset<0 || $offset >= length($seq));	
            # the situation sometimes should be change dependent on Flag of different alignmentor
			if($read_strand & 0x10)
            {
				my $qual_score = ord(substr($qual_string,$offset+1,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score);
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=$rcTable{substr($seq,$offset+1,1)};			
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;	
			}else{
				my $qual_score = ord(substr($qual_string,$offset,1))-$phred;
				next if($readHapInfo{$read_id}->{$CpG_position}->{"qual"} && $readHapInfo{$read_id}->{$CpG_position}->{"qual"} > $qual_score);
				$readHapInfo{$read_id}->{$CpG_position}->{"base"}=substr($seq,$offset,1);
				$readHapInfo{$read_id}->{$CpG_position}->{"qual"}=$qual_score;

			}	
			#print "$offset\t", $readHapInfo{$read_id}->{$CpG_position}, "\n";
		}	
		$read_start_pos_table{$read_id} = $read_start if(!$read_start_pos_table{$read_id} || $read_start_pos_table{$read_id}<$read_start);
	}		
	
	#Use read start as the UMI to collapse multiple clonal reads.
	my @read_list = keys(%readHapInfo);		
	my %unique_read_base_info;
	foreach my $read_id (@read_list){		
		my $UMI = $read_start_pos_table{$read_id};
		foreach my $CpG_position (keys(%{$readHapInfo{$read_id}})){
			next if(!$readHapInfo{$read_id}->{$CpG_position}->{"base"});
			#print "$CpG_position\n";
			#print $readHapInfo{$read_id}->{$CpG_position}->{"base"}, "\n";
			#print $readHapInfo{$read_id}->{$CpG_position}->{"qual"}, "\n";
			$unique_read_base_info{$UMI}->{$CpG_position}->{$readHapInfo{$read_id}->{$CpG_position}->{"base"}}+=$readHapInfo{$read_id}->{$CpG_position}->{"qual"};
		}
	}
	
	#Derive the consensus haplotype string based on multiple clonal reads.
	foreach my $UMI (keys(%unique_read_base_info)){
		my $hap_string="";
		foreach my $CpG_position(@CpG_positions){
			my $base = "N";
			my $best_qual = 0;
			if($unique_read_base_info{$UMI}->{$CpG_position}){
				foreach my $base_call (keys(%{$unique_read_base_info{$UMI}->{$CpG_position}})){
					next if($unique_read_base_info{$UMI}->{$CpG_position}->{$base_call} <= $best_qual);
					$best_qual = $unique_read_base_info{$UMI}->{$CpG_position}->{$base_call};
					$base = $base_call;
				}
			}					
			$hap_string=$hap_string.$base;
		}	
		$hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string}++ if($hap_string);
	}

	#report haplotype strings only on the positions with valid calls.
	foreach my $hap_string(keys %{$hapInfoTable{$target_id}->{"hap_counts"}}){
		my @valid_positions;
		my $valid_hap="";
		for(my $i=0; $i<length($hap_string); $i++){
			my $allele = substr($hap_string,$i,1);
			next if($allele eq "N");
			$valid_hap = $valid_hap.$allele;
			push(@valid_positions,$CpG_positions[$i]);
		}
		next if(scalar(@valid_positions)<1);
		next if($valid_hap =~ /[AG]/);
		print $target_id,"\t$valid_hap\t", $hapInfoTable{$target_id}->{"hap_counts"}->{$hap_string},"\t",join(",",@valid_positions), "\n";
	}
}
}
	
