 #!/usr/bin/perl

 # creat whole encode TFBS file to seperate files by TF

 use strict;
 use Cwd;
 my $dir=getcwd;
 chdir $dir or die "can't change into $dir\n";
 my $input="wgEncodeRegTfbsClusteredWithCellsV3.bed";

 open F,"$input";
 while(<F>){
 my ($chr,$start,$end,$tf,undef)=split /\t/;
 open OUT,">>encode.$tf.hg19.bed";
 print OUT $_;
 }
