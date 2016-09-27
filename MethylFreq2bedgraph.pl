# This script could transfer the MethylFreq files from Dinh to bedgraph files so that I can upload them to UCSC. 

#!/usr/bin/perl -w

use strict;

my $input = $ARGV[0];
my($sampleID,undef) = split /\./,$input;
my $minDepth=$ARGV[1];

print "USAGE: perl $0 MethylFreq 5 > prefix.bedgraph\n";

$minDepth = 5 if(!$minDepth);

print "track type=bedgraph name=\"$sampleID\" visibility=full color=20,150,20 altColor=150,20,20 windowingFunction=mean\n";

my %methylTable;
sub main(){
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
        report_methylFreqBED();
}

sub report_methylFreqBED(){
        my $cur_chr = "NA";
        foreach my $chr(sort keys(%methylTable)){
                foreach my $pos(sort {$a<=>$b} keys %{$methylTable{$chr}}){
                        next if($methylTable{$chr}->{$pos}->{'CT'}<$minDepth);
                        my $methylLevel = sprintf("%4.3f", $methylTable{$chr}->{$pos}->{'C'}/$methylTable{$chr}->{$pos}->{'CT'});
						my $start=$pos;
						my $end=$pos+1;
                        print "$chr\t$start\t$end\t$methylLevel\n";
                }
        }
}

main();
