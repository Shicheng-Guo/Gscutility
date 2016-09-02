#!/usr/bin/perl
# search the postions for the CpGs in human genome

my $fa=shift @ARGV;
open F,$fa;
my $seq;
while(<F>){
chomp;
$seq .=$_;
}
close F;

sub match_positions {
    my ($regex, $string) = @_;
    return if not $string =~ /($regex)/;
    return (pos($string), pos($string) + length $1);
}
sub all_match_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /($regex)/g) {
        push @ret, [(pos($string)-length $1),pos($string)-1];
    }
    return @ret
}

my $regex='CG';
my $string=$seq;
my $cgap=3;    
my @pos=all_match_positions($regex,$string);
my @hgcg;
foreach my $pos(@pos){
    push @hgcg,@$pos[1];
}
foreach my $i(0..($#hgcg)){
my $pos=$hgcg[$i]-1;  # transfer to 0-based coordinate 
print "$chr\t$pos\n"; 
}


sub USAGE(){
print "
------------------------------------------
perl $0 chr12.fa 

Find all the CpG postion (0-base) from genome Fa files (Download from UCSC)


------------------------------------------
"
}

