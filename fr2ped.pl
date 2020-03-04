use strict;
use Cwd;
chdir getcwd;
open F,shift @ARGV;
my $i=1;
while(<F>){
next if !/ZS/;
my($snp,$sam,$rs,$gc,$chr,$pos,$a1,$a2,undef)=split/\s+/;
if($i eq 1){
print "$sam $sam 0 0 0 0 $a1 $a2";
}else{
print " $a1 $a2";
}
$i++;
}
print "\n";
