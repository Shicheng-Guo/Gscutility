#!/usr/bin/perl

package readDataset;
=head1 readDataset.pl

this script is an example of reading the samples correlation dataset, 
it has been written with clarity in mind, rather than efficiency

it requires the probes list as index file
=cut

use strict;
use warnings;

=head2 readIndex(filepath)

Read an input text file, returns a hash with each line as hash keys 
and hash values as index number ( 0 based ), this hash will serve
as index in most of our code.
We assume that the input file contains no duplicate, and would fail 
if the input contain such duplicate

=cut
sub readIndex($) {
	my ( $input ) = @_;
	my $fh;
	my %result;
	my $i=0;
	open( $fh, "<", $input ) or die( "error reading [$input] : $!" );
	while(<$fh>) {
		my $item = $_;
		chomp($item);
		if( defined( $result{$item} ) ) {
			die( "duplicate entry in the input [$input] : [$item]" );
		}
		$result{$item} = $i++;
	}
	close( $fh );
	return(%result);
}

=head2 readDatasetLinear(filepath,indexHash)

Read all the dataset and for each entry print :

	probe1	probe2	correlation

=cut
sub readDatasetLinear($%) {
	my ( $input, %index ) = @_;

	my @keys=keys( %index );
	#my $count=scalar( @keys );

	my $fh;
	open( $fh, "<", $input ) or die( "error reading [$input] : $!" );

	my $buffer='';
	my $read=0;
	while(!eof($fh)) {
		my $offset=tell($fh);

		$buffer='';
		$read=read($fh,$buffer,2,0);
		die("error reading [$input]@$offset") unless( $read == 2 );
		my $key1 = unpack('s',$buffer);

		$buffer='';
		$read=read($fh,$buffer,2,0);
		die("error reading [$input]@$offset") unless( $read == 2 );
		my $key2 = unpack('s',$buffer);

		$buffer='';
		$read=read($fh,$buffer,4,0);
		die("error reading [$input]@$offset") unless( $read == 4 );
		my $score = unpack( 'f', $buffer );

		my $probe1 = $keys[$key1];
		my $probe2 = $keys[$key2];
		print( "$probe1\t$probe2\t$score\n" );
	}
	close( $fh );
}

=head2 readDatasetPair(filepath,probe1,probe2,indexHash) 

Read a specific score from the input file corresponding to the 
query pair of probe.

=cut
sub readDatasetPair($$$%) {
	my ( $input, $probe1, $probe2, %index ) = @_;

	if( $probe1 eq $probe2 ) {
		return( 1.0 );
	}

	my $key1 = $index{$probe1};
	my $key2 = $index{$probe2};

	die("error : invalid probe specified ") unless ( defined( $key1 ) and defined( $key2 ) ); 

	if($key1 > $key2) {
		($key1,$key2) = ($key2,$key1)
	}

	my @keys = keys( %index );

	my $n = scalar( @keys );

	my $offset = 0;

	# get linear index
	my $k = ( $n * ( $n - 1) / 2 ) - ( ( $n - $key1 ) * ( $n - $key1 - 1) / 2 ) + $key2 - $key1 - 1 ;
	my $fh;
	open( $fh, "<", $input ) or die( "error reading [$input] : $!" );
	my $score =  0.0;
	if( seek( $fh, $k * 8 + 4, 0) == 1 ) {
		my $buffer='';
		my $read = read( $fh, $buffer, 4 , 0);
		die( "error reading [$input] : [$read] : [$!]" ) unless( $read == 4 );
		$score = unpack( 'f', $buffer );
	}else{
		die("error reading [$input] : $!");
	}

	close( $fh );
	return $score;
}

sub showUsage() {
	print <<"END";
Usage : 
	perl readDataset.pl \$indexFile \$inputFile [\$probe1 \$probe2]

Will print the dataset in the input file based on the index file.
if a pair of probes are specified, the program will print 
only the score of that probe pair.

END
}

#####################################################################

my ( $indexFile, $inputFile, $probe1, $probe2 ) = @ARGV;

unless ( defined( $indexFile ) and defined($inputFile) ) {
	showUsage();
	die();
}

my %index=readIndex( $indexFile );

if( defined( $probe1 ) and defined( $probe1 ) ) {
	my $score=readDatasetPair($inputFile,$probe1,$probe2,%index);
	print( "$probe1\t$probe2\t$score\n" );
}else{
	readDatasetLinear($inputFile,%index);
}

__END__
