################################################################################
## the script used to merge all Rockhopper_Results/NC_007793_transcripts.txt
## project: RNA-seq data to reveal novel response mechanism to bacterial within host wound tissues
## Shicheng Guo, Shihcheng.Guo@Gmail.com
## 02/11/2020

use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;
my $chr=shift @ARGV;
my @file=glob("*/NC_007793_transcripts.txt");

my %data;
my %key;
my %sam;
my %name;
foreach my $sam(@file){
my $input="/home/guosa/hpc/project/RnaseqBacterial/extdata/rnaseq/Rockhopper_Results/".$sam;
my ($sam,$name)=split /\//,$sam;
$sam{$sam}=$sam;
#print "$input\n";
open F,$input || die "$input cannot be open!\n";
while(<F>){
        next if /RelativeActivity/;
        next if /Transcription/;
        my $value;
        chomp;
        my @line=split /\t/;
        #my $check=join("+",@line);
        #print "$check\n";
        next if $line[6]=~/predicted/;
        $name{$line[6]}="$line[5]\t$line[7]";
        next if $line[6] eq "";
        my $key="$line[6]";
        $value=$line[$#line];
        $key{$key}=$key;
        $data{$key}{$sam}=$value;
        }
}
my @sam;
foreach my $sam(sort keys %sam){
        push(@sam,$sam);
}
my $head=join("\t",@sam);
print "\t$head\n";

foreach my $key(sort keys %data){
        print "$key";
        foreach my $sam (sort keys %sam){
        if( ! exists $data{$key}{$sam}||! defined $data{$key}{$sam}||$data{$key}{$sam}=~/NA/){
        print "\t0";
        }else{
        my $R=sprintf("%.2f",$data{$key}{$sam});
        print "\t$R";
        }
        }
        print "\t$name{$key}\n";
}
