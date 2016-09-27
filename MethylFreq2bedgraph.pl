# The script will transfer all the MethylFreq (Dinh) files to bedGraph in current fold

#!/usr/bin/perl -w
use strict;
use Cwd;

my $dir=getcwd;
chdir $dir;

my @file=glob("*.methylFreq");

print "USAGE: perl $0 \n";
print "The script will transfer all the MethylFreq files to bedGraph in current fold\n";

foreach my $input (@file){ 
	my $minDepth=1;
	my($sampleID,undef) = split /\./,$input;
	my $header= "track type=bedGraph name=\"$sampleID\" visibility=full color=20,150,20 altColor=150,20,20 windowingFunction=mean\n";
	my $output="$sampleID.bedGraph";
	my %methylTable;
	open F,$input;
        while(my $line = <F>){
                chomp($line);
                my @fields = split(/\s+/, $line);
                my $strand = $fields[2] eq 'W' ? '+' : '-';
                my %alleleCounts;
                my $CT_counts;
                for(my $i=5; $i<scalar(@fields); $i+=2){
                        $alleleCounts{$fields[$i]}=$fields[$i+1];
                        $CT_counts += $fields[$i+1] if($fields[$i]=~ /[CT]/);
                }
                next if(!$CT_counts || $CT_counts/$fields[3] < 0.9);
                my $index=$fields[0] . ":" . $fields[1];
                $alleleCounts{'C'} =0 if(!$alleleCounts{'C'});
                $methylTable{$fields[0]}->{$fields[1]}->{'C'} +=  $alleleCounts{'C'} ;
                $methylTable{$fields[0]}->{$fields[1]}->{'CT'} += $CT_counts;
        }
        &report_methylFreqBED(\%methylTable,$minDepth,$header,$output);
}

sub report_methylFreqBED(){
        my $methylTable=shift @_;
        my $minDepth=shift @_;
        my $header=shift @_;
	my $output=shift @_;
	my %methylTable=%{$methylTable};
	open OUT,">$output";
        my $cur_chr = "NA";
	print OUT $header;
        foreach my $chr(sort keys(%methylTable)){
                foreach my $pos(sort {$a<=>$b} keys %{$methylTable{$chr}}){
                        next if($methylTable{$chr}->{$pos}->{'CT'}<$minDepth);
                        my $methylLevel = sprintf("%4.3f", $methylTable{$chr}->{$pos}->{'C'}/$methylTable{$chr}->{$pos}->{'CT'});
						my $start=$pos;
						my $end=$pos+1;
                        print OUT "$chr\t$start\t$end\t$methylLevel\n";
                }
        }
	close OUT;
}

