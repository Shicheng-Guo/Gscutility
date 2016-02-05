#!/usr/bin/perl 

#Tyler K. Chafin; 29-Feb-2016 
#tkchafin@uark.edu 

use strict; 
use warnings; 
use Getopt::Std;

my $usage = "
This script functions to convert a 2-column snp file to input for Arlequin. Written by Tyler K. Chafin for Bifang He as per request on ResearchGate forum. 

Author: Tyler K. Chafin - tkchafin\@uark.edu
Last Modified: 22-Dec-15

Usage: $0 </path/to/.snp> <#cols to skip> 

"; 


#Read command line arguments
my %opts;
getopts( 'f:hi:c:i:p', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
}

my ($file, $cols, $id, $pop) = &parseArgs(\%opts);
#print "$file\n$cols\n$id\n";

print "\nOperating on file: $file\n"; 
print "Skipping columns before SNP matrix: $cols\n"; 
print "Looking for ID in column: $id\n"; 

my @line; 
my %contents; 
my $count = 0; 
open (INPUT, $file) || die "Cannot open $file: $!\n\n"; 
while (<INPUT>){ 
  chomp;   
  my $temp; 
  my @linesnps; 
  if ($_ =~ /^\s*$/){
    next; 	  
  }else{
    @line = split /\s+/, $_; 
    if ($id > 0){ 
      $temp = $line[$id-1]; 
    }else{
      $temp = $count;
    }
    #Subset array
    if ($contents{$temp}){
      print "Warning: Sample id $temp occurs more than once! Are you sure you defined the correct column for the -i option?\n\n"; 
    }
    $contents{$temp} = [@line[$cols..$#line]];

    $count++; 
  }   
}
close INPUT;

open (OUT, ">$file.arp") || die "Cannot open $file.arp: $!\n";
print "Writing output file:$file.arp\n";
print OUT "[PROFILE]\n"; 
print OUT "  DataType = DNA\n";
print OUT "  NbSamples = 1\n";
print OUT "  GameticPhase=0\n";
print OUT " LocusSeparator=NONE\n";
print OUT "  GenotypicData = 1\n";
print OUT "\n[DATA]\n";
print OUT "[[Samples]]\n";

print OUT "  SampleName=\"1\"\n";
print OUT "  SampleSize=$count\n";
print OUT "  SampleData={\n";

my $line1;
my $line2;
foreach my $key ( sort { $a<=>$b } keys %contents){
  $line1 = "    ";
  $line2 = "    ";
  my $i= 0;
  $line1 .= "$key 1\t" . join("", grep {$i++ % 2}  @{$contents{$key}}); $i=1;
  $line2 .= "\t" . join("", grep {$i++ % 2} @{$contents{$key}}); 
  print OUT "$line1\n$line2\n";
}
print OUT"  }\n";
close OUT; 
print "Done!\n\n";
exit; 

############################## SUBROUTINES #############################
sub help{
  print"
This script functions to convert a 2-column snp file to input for Arlequin. Written by Tyler K. Chafin for Bifang He as per request on ResearchGate forum.

Author: Tyler K. Chafin - tkchafin\@uark.edu
Last Modified: 22-Dec-15

Usage: $0 -f <input file> -i <ID column> -c <columns before snp>

Options 
	-f	- Input file (white-space separated SNP genotypes; 2 column each)
	-i	- Which column contains individual IDs 
		  Default = 0; meaning no identifier is given. 
	-c	- How many columns to skip
		  Defauls = 0 columns will be skipped
	-h	- Kills program and prints this help message 

";
exit; 
}

sub parseArgs{

  my ($params) = @_; 
  my %opts = %$params;

  # set default values for command line arguments
  my $file = $opts{f} or &help;
  my $cols = $opts{c} or $cols = 0;
  my $id = $opts{i} or $id = 0;

  return( $file, $cols, $id );

}
