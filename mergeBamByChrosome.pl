#!/usr/bin/perl
# merge bam file by different chrosome(1-22, M, X, Y)
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-01-19

# Bam File: /media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles
# OUT Dir: /home/shg047/oasis/Estellar2016/bam

use strict;
use Cwd;
my $bamdir="/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles";
chdir $bamdir;
my @file=glob("*bam");
my $outdir="/home/shg047/oasis/Estellar2016/bam";
my %sam;
foreach my $file(@file){
        my ($sam,undef)=split /\./,$file;
        $sam{$sam}=$sam;
}
foreach my $sam(sort keys %sam){
	open OUT, ">$outdir/$sam.bam.merge.sh";
	print OUT "cd $bamdir\n";
	my @file=glob("$sam*bam");
	my $file=join(" ",@file);
	print "$sam.bam.merge.sh\n";
    print OUT "samtools view -H $file[0]>$outdir/$sam.header\n";
    print OUT "samtools cat -h $outdir/$sam.header -o $outdir/$sam.bam $file\n";
	print OUT "samtools index $outdir/$sam.bam\n";
	close OUT ;			
	system("sh $outdir/$sam.bam.merge.sh &");	
	#system("sh $outdir/$sam.bam.merge.sh");	
}










