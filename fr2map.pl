use strict;
use Cwd;
chdir getcwd;
open F,shift @ARGV;
while(<F>){
next if !/ZS/;
my($snp,$sam,$rs,$gc,$chr,$pos,$a1,$a2,undef)=split/\s+/;
print "$chr $rs 0 $pos\n";
}

