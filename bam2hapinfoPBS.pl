#!/usr/bin/perl -w
use strict;
use Cwd;

my $bam2hapinfo = "/home/shg047/bin/mergedBam2hapInfo.pl";
my $sample_info_file = $ARGV[0];
my $submit=$ARGV[1];

my %sample_group;
die &USAGE if scalar(@ARGV)<1;

open(INFILE, "$sample_info_file")||die("Error in opening file $sample_info_file.\n");
while(my $line = <INFILE>){
                chop($line);
                my ($id,$bam_file, $target_bed) = split(/\t/, $line);
                next if(!$id || !$bam_file);
                $sample_group{$id}->{"bam_file"}=$bam_file;
                $sample_group{$id}->{"target_bed"}=$target_bed;
}
close(INFILE);

foreach my $id (sort keys(%sample_group)){
    next if $id=~/LambdaNEB/i;
    my $job_file_name = $id . ".job";
    my $status_file = $id.".status";
        my $hapInfo_file = $id.".hapInfo.txt";
        my $curr_dir = `pwd`;
        open(JOB_FILE, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
        print JOB_FILE "#!/bin/csh\n";
        print JOB_FILE "#PBS -q glean\n";
        print JOB_FILE "#PBS -l nodes=1:ppn=1\n";
        print JOB_FILE "#PBS -l walltime=8:00:00\n";
        print JOB_FILE "#PBS -o ".$id.".log\n";
        print JOB_FILE "#PBS -e ".$id.".err\n";
        print JOB_FILE "#PBS -V\n";
        print JOB_FILE "#PBS -M shihcheng.guo\@gmail.com \n";
        print JOB_FILE "#PBS -m abe\n";
        print JOB_FILE "#PBS -A k4zhang-group\n";
        print JOB_FILE "cd $curr_dir\n";
        my $cmd = "$bam2hapinfo ".$sample_group{$id}->{"target_bed"}." ".$sample_group{$id}->{"bam_file"}." > $hapInfo_file";
        print JOB_FILE "$cmd\n";
        close(JOB_FILE);
        print "Job file $job_file_name created.\n";

        system("qsub $job_file_name");
}

sub USAGE{
	print "USAGE: perl $0 sample_info_file submit\n";
	
	print "Run to creat SamInfo: perl samInfoPrep4Bam2Hapinfo.pl ./\n";
	print '
Sample_Info_File Format: 
ENCFF000LUN_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUN_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUP_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUP_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUQ_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUQ_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUT_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUT_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUU_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUU_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LUV_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LUV_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVA_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVA_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVB_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVB_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVE_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVE_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVF_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVF_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVI_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVI_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
ENCFF000LVJ_trimmed     /home/shg047/oasis/Haib/bam/ENCFF000LVJ_trimmed.fq.gz_bismark_bt2.sort.bam      /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
'
}
