
#!/usr/bin/perl -w

# Transfer Bam to Fastq with samtools command
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-02-19

use strict;
use Cwd;

chdir getcwd;
my %sam;
my @fastq=glob("*fastq.gz");
foreach my $fastq(@fastq){
        my($sample,undef)=split /[_.]/,$fastq;
        $sam{$sample}=$sample;
}

foreach my $sam (sort keys %sam){
        print "$sam\_1.fastq.gz\t$sam\_2.fastq.gz\n"
}
