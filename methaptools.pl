#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
chdir getcwd;
use Getopt::Long;
use File::Basename;
use IO::Handle;
$|++;

# Extract Pileup within specfic genomic regions (CDS, Exon)
# Set Pileup and GenomicInterval regions
# Contact: Shicheng Guo
# Version 1.3
# Update: 2017-04-5

#chdir "C:\\Users\\shicheng\\Downloads\\hapinfo2mhl";
my $haptools="v0.16.3_dev";
my ($input,$output,$bed,$help)=commandopt();
#my $input="Indx04.sortc.bam.hapInfo.txt.head";
#my $output="Indx04.hap";
#my $bed="input.bed";


# step 1. read bed regions from bed files
open F,$bed;
my $window=5000;
my $bedthreshold=50000;
my %hash;
my %bed;
my %HapinfoRegion;
my %hapCountMmatrix;

while(<F>){
    chomp;
    next if /^\s+$/;
    my ($chr,$start,$end,$id,$gene,$block)=split/\s+/;
	my $len=$end-$start;
    if($len>$bedthreshold){
	warn "Genomic regions from $bed should < $bedthreshold! please remove large regions!";  
	next;
	}
    my $bin1=int($start/$window);
    my $bin2=int($end/$window);
    my $target="$chr:$start-$end";
    $bed{$target}=$target;
    foreach my $i($bin1..$bin2){
        $hash{$chr}{$i-1}{$target}=$target;
        $hash{$chr}{$i}{$target}=$target;
        $hash{$chr}{$i+1}{$target}=$target;
	}
}
close F;

# step 2. pop haptype into memory
open F1,$input || die "cannot open $input\n";
print "\nStart to reading haplotype files....\n";
my($filename, $dirs, $suffix) = fileparse($input);
while(<F1>){
chomp;
next if /^\s*$/;
#print "$_\n";                           # bug
my ($probeID,$hapString,$Hapcount,$cpgPos)=split /\s+/;
next if length($hapString)<1;            # option
my ($chr,undef,undef)=split/[:-]/,$probeID;
my @tmp=split/,/,$cpgPos;
my $bin=int($tmp[0]/$window);

my $readstart=$tmp[0];
my $readend=$tmp[$#tmp];
if(exists $hash{$chr}{$bin}){
	foreach my $bed(sort keys %{$hash{$chr}{$bin}}){
		# print "\t$chr\t$bin\t$bed\t";   # debug
		my ($chr,$start,$end)=split/[:-]/,$bed;
		my $HapString;                    # keep shorted haplotypes from substr
		if($readstart>=$start){
		my @pos;
		foreach my $i(0..$#tmp){
		push @pos,$i if ($tmp[$i]<=$end);
		# print "$i," if ($tmp[$i]<=$end and $tmp[$i]>=$start);	
		}
		if(scalar(@pos)==1){
		$HapString=substr($hapString,$pos[0],1);
		$cpgPos=substr($cpgPos,$pos[0],1);
		$hapCountMmatrix{$bed}->{$cpgPos}->{$HapString}+=$Hapcount;
		}elsif(scalar(@pos)>1){
		# my $len=scalar(@pos);               # debug
		# print "$len::$hapString\t";         # debug
		$HapString=substr($hapString,$pos[0],scalar(@pos));
		$hapCountMmatrix{$bed}->{$cpgPos}->{$HapString}+=$Hapcount;
		}
		# my $len=scalar(@pos);               # debug
		# print "$len::$hapString\t";         # debug
		}else{
			my @pos;
			foreach my $i(0..$#tmp){
			push @pos,$i if ($tmp[$i]<=$end and $tmp[$i]>=$start);	
			# print "$i," if ($tmp[$i]<=$end and $tmp[$i]>=$start);	
			}
			if(scalar(@pos)==1){
			my $len=scalar(@pos);                          # debug
			# print "$_\t$bed\t$hapString\t$pos[0]\t$len\t";
			$HapString=substr($hapString,$pos[0],1);
			# print "substr:$HapString\n";
			$hapCountMmatrix{$bed}->{$cpgPos}->{$HapString}+=$Hapcount;
			}elsif(scalar(@pos)>1){
			# my $len=scalar(@pos);             # debug
		    #print "$len::$hapString\t";       # debug
			#print "$_\t$bed\t$hapString\t$pos[0\t$len\n";
			$HapString=substr($hapString,$pos[0],scalar(@pos));
			$hapCountMmatrix{$bed}->{$cpgPos}->{$HapString}+=$Hapcount;
			}
			# my $len=scalar(@pos);             # debug
			# print "$len::$hapString\t";       # debug
		}		
	}
	# print "\n";
}
}
close F1;
print "Haplotype information: $input loading completed, caculating is on the way!";
open OUT1,">$output.hap";
print OUT1 "Postion\tMF\tMHL\tCCT\tHom\tNCC\tNTT\tNCT\n";
my $MF=&hapinfo2mf(\%hapCountMmatrix);
my ($TEMP)=&hapinfosummary(\%hapCountMmatrix);
my ($MHL,$CCP)=@{$TEMP};
foreach my $bed(sort keys %{$MHL}){
	print OUT1 "$bed\t$MF->{$bed}\t";
	print OUT1 "$MHL->{$bed}\t";
	foreach my $ccp(@{$CCP->{$bed}}){
	print OUT1 "$ccp\t";
	}
	print OUT1 "\n"
}
close OUT1;

# Step 3. Summary and Calculations(MHL,CCP,CNT)
sub hapinfo2mf(\%){
	# $hapCountMmatrix{$bed}->{$cpgPos}->{$hapString}
	my %MF;
	my ($hapCountMmatrix,undef)=@_;
	my %hapCountMmatrix=%{$hapCountMmatrix};
	foreach my $probeID(sort keys %hapCountMmatrix){
		my $HapStrings;
		my ($C,$T,$mf);
		foreach my $cpgPos(sort keys %{$hapCountMmatrix{$probeID}}){
			foreach my $hapString(sort keys %{$hapCountMmatrix{$probeID}{$cpgPos}}){
				next if($hapString =~ /[NAG]/i);
				foreach (my $i=1;$i<=$hapCountMmatrix{$probeID}->{$cpgPos}->{$hapString};$i++){
				$HapStrings .=$hapString;
				}
			}
		}
		my $C_number = () = $HapStrings =~ /C/gi;
		my $T_number = () = $HapStrings =~ /T/gi;
        $mf = ($C_number/($C_number+$T_number));
	    # print "$probeID\tC:$C_number\tT:$T_number\tmf:$mf\t$HapStrings\n";
	    $MF{$probeID}=$mf;	
		}
	return \%MF;
}
	
sub hapinfosummary(\%){
my ($hapCountMmatrix,undef)=@_;
my %hapCountMmatrix=%{$hapCountMmatrix};
my (%mhl_load_matrix, %HapCCP,@output);
foreach my $probeID (sort keys(%hapCountMmatrix)){
	my %k_mer_counts;
	my $mc_hap_load=0;
	foreach my $cpgPos(sort keys %{$hapCountMmatrix{$probeID}}){
		foreach my $hapString (keys(%{$hapCountMmatrix{$probeID}{$cpgPos}})){
			for(my $word_size = 1; $word_size<=length($hapString); $word_size++){
				next if($word_size>10);
				for(my $i=0; $i<=length($hapString)-$word_size; $i++){
					my $sub_hapString = substr($hapString,$i,$word_size);
					next if($sub_hapString =~ /[NAG]/i);
					$k_mer_counts{$word_size}->{$sub_hapString}+=$hapCountMmatrix{$probeID}->{$cpgPos}->{$hapString};
				}
			}
		}
	}
	# my $norm_factor=0;
	# print "$probeID\t"; # debug
    my $mhl=&kmerc2mhl(\%k_mer_counts);
    my $ccp=&kmerc2ccp(\%k_mer_counts);
	$mhl_load_matrix{$probeID}=$$mhl;
	push @{$HapCCP{$probeID}},@{$ccp};
	
}	
    @output=(\%mhl_load_matrix,\%HapCCP);
	return(\@output);
}

sub kmerc2mhl(\%){
my ($k_mer_counts,undef)=@_;
my %k_mer_counts=%{$k_mer_counts};
my $mc_hap_load;
my $norm_factor=0;
my %mhl_load_matrix;
my @unmethylated_haps= ("T"x1,"T"x2,"T"x3,"T"x4,"T"x5,"T"x6,"T"x7,"T"x8,"T"x9,"T"x10,"T"x11,"T"x12,"T"x13,"T"x14);
my @methylated_haps  = ("C"x1,"C"x2,"C"x3,"C"x4,"C"x5,"C"x6,"C"x7,"C"x8,"C"x9,"C"x10,"C"x11,"C"x12,"C"x13,"C"x14);
foreach my $word_size (sort keys(%k_mer_counts)){
	$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$unmethylated_haps[$word_size-1]});
	$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}=0 if(!$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]});
	my $total_count=0;
	foreach my $allele (keys(%{$k_mer_counts{$word_size}})){
		$total_count+=$k_mer_counts{$word_size}->{$allele};
		#print "$sample_name\t$allele\t$k_mer_counts{$word_size}->{$allele}\n";# debug
	}
	next if($total_count<1);
	my $mh_fraction = $k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}/$total_count;
	my $weight = $word_size;
    #print "$sample_name\t$word_size\t$methylated_haps[$word_size-1]\t$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}\t$total_count\t$weight\t$mh_fraction\n";# debug
	$mc_hap_load += $weight*$mh_fraction;
	#print "$k_mer_counts{$word_size}->{$methylated_haps[$word_size-1]}/$total_count x $weight\t";# debug
	$norm_factor+=$weight;		
	}
	next if(!$norm_factor);
	$mc_hap_load/=$norm_factor;
	#print "=$mc_hap_load/$norm_factor=$mc_hap_load\n"; # debug
	$mc_hap_load=sprintf("%.4f",$mc_hap_load);
	return(\$mc_hap_load)
}

sub kmerc2ccp(\%){
# $k_mer_counts{$word_size}->{$sub_hapString}+=$hapCountMmatrix{$probeID}->{$cpgPos}->{$hapString};
my ($k_mer_counts,undef)=@_;
my %k_mer_counts=%{$k_mer_counts};
my (%HapCCP,%HapinfoCCP,$C,$T,$N,$ratio1,$ratio2);
foreach my $word_size (sort keys(%k_mer_counts)){
	foreach my $hapString(sort keys %{$k_mer_counts{$word_size}}){
		if($hapString eq "C"x $word_size){
    		$HapinfoCCP{'C'}+=$k_mer_counts{$word_size}->{$hapString};
    	}elsif($hapString eq "T"x $word_size){
    		$HapinfoCCP{'T'}+=$k_mer_counts{$word_size}->{$hapString};
    	}else{
    		$HapinfoCCP{'N'}+=$k_mer_counts{$word_size}->{$hapString};
		}
	}
}
$C= defined $HapinfoCCP{'C'} ? $HapinfoCCP{'C'} : "0"; 
$T= defined $HapinfoCCP{'T'} ? $HapinfoCCP{'T'} : "0"; 
$N= defined $HapinfoCCP{'N'} ? $HapinfoCCP{'N'} : "0"; 
my $total=$C+$T;
if($total == 0){
$ratio1="NA";
$ratio2="NA";
}else{
$ratio1=sprintf("%.3f",$C/($C+$T));
$ratio2=sprintf("%.3f",1-($N)/($C+$T+$N));
}
my @result=($ratio1,$ratio2,$C,$T,$N);
return(\@result);
}

	
sub commandopt{
	my ($help,$input,$output,$bed);
	my $command_line=GetOptions ( "input=s"=> \$input,"output=s"=> \$output,"bed=s"=> \$bed,"help"=>\$help);
	return($input,$output,$bed,$help);
}

sub USAGE{
print "Usage: perl methhaptools.pl --input SRX209456.hapInfo.txt --bed interest.bed --output test\n";
print "Methylation Haplotype Related Ratio and Counts in for Each Haplotype Files\n";

print '
Format For: Interest genomic interval file: 
chr10:76532564-76532591		CTCC	19	76532564,76532564,76532564,76532564
chr10:76532564-76532591		CTTC	19	76532564,76532564,76532564,76532564

Format For: Interest genomic interval file: 
chr10   76532564        76532591        Col18a1 Exon
chr10   76532243        76532564        Col18a1 Intron
';
}





