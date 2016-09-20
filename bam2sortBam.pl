
#!/usr/bin/perl -w

# Bam to sort Bam with samtools in TSCC
# Run the script to the Bam directory
# Contact: Shicheng Guo
# Version 1.3
# Update: 2016-01-19

use strict;
use Cwd;
die USAGE() if scalar(@ARGV)<1;
my $dir=getcwd;
my $sort=shift @ARGV;
my $submit=shift @ARGV;

my @file=glob("$dir/*.bam");
mkdir "../sortBam" if ! -e "../sortBam";
my $ppn=6;
my $walltime="168:00:00";
my $queue="hotel"; # hotel

foreach my $file (@file){
    my @tmp=split /\//,$file;
    my ($sample,undef)=split /.bam/,$tmp[$#tmp];
    my $job_file_name = $sample.".job";
    my $status_file = $sample.".status";
    my $hapInfo_file = $sample.".hapInfo.txt";
    my $curr_dir = $dir;
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $job_file_name\n";
    print OUT "#PBS -q $queue\n";
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    print OUT "#PBS -l walltime=$walltime\n";
    print OUT "#PBS -o ".$sample.".log\n";
    print OUT "#PBS -e ".$sample.".err\n";
    print OUT "#PBS -V\n";
    print OUT "#PBS -M shicheng.guo\@gmail.com \n";
    print OUT "#PBS -m abe\n";
    print OUT "#PBS -A k4zhang-group\n";
    print OUT "cd $curr_dir\n";
    if($sort eq 'sortc'){
    print OUT "samtools sort -@ 6 -o ../sortBam/$sample.sortc.bam $file\n";
    print OUT "samtools index ../sortBam/$sample.sortc.bam\n";
    }elsif($sort eq 'sortn'){
    print OUT "samtools sort -n -@ 6 -o ../sortBam/$sample.sortn.bam $file\n";
    print OUT "samtools index ../sortBam/$sample.sortn.bam\n";
    }
    if($submit eq 'submit'){
    system("qsub $job_file_name");
    }
}

sub USAGE{
        print "\nUSAGE: perl $0 sortc|sortn submit\n\n";
        print "Example: perl $0 /media/LTS/shg047/bam sortc submit\n";
}
