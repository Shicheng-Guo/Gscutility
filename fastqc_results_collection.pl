use strict;
use warnings;

my $usage="$0 <listf> <outf>\n
<listf>  the fastQC output summary file list in the format:
\t\tSampleID<tb>InsertSize<tb>fastqc1<tb>fastqc2, the InsertSize and fastqc2 is omitted for SE
<outf>   the statistics result file


\n";

die $usage if @ARGV<2;

my($listf,$outf)=@ARGV;

my $wc = `cat $listf |head -n1 |sed 's/\\t/\\n/g'|wc -l`;chomp $wc;
my $ps;
if($wc > 2){
        $ps = "PE";
}else{
        $ps = "SE";
}

my(@info,%res);
open I,$listf;
open O,">",$outf;
if($ps eq "PE"){
        print O "Sample\tInsertSize\tReadLength\tGC%\tRawReads\tRawBases\tQ20%\tQ30%\n";
}else{
        print O "Sample\tReadLength\tGC%\tRawReads\tRawBases\tQ20%\tQ30%\n";
}
while(<I>){
        my ($total,$ratio,%posNum,$readLen,%res,$readsNum1,$readsLen1,$fastqcRes);
        chomp;
        $fastqcRes=$_;
        @info=split /\s+/;
        $res{InsertSize}=$info[1];
        my $samp=$info[0];
        if($ps eq "PE"){
                open I2,$info[2];
                open I3,$info[3];
        }else{
                open I2,$info[1];
        }
        while(<I2>){
                chomp;
                if(/Total\sSequences\s(\d+)/){
                        $total=$1;
                        $posNum{20}=$total;
                        $posNum{30}=$total;
                        $res{RawReads}= $total;
                        $readsNum1=$1;
                }
                if(/Sequence\slength\s(\d+)[-_](\d+)/){
                        warn "can not correctly statistic the base number, because the reads length is not the same, use the shortest one\n$_\n";
                }
                if(/Sequence\slength\s(\d+)/){
                        $res{ReadLength}=$1."_".$1 if $ps eq "PE";
                        $res{ReadLength}=$1 if $ps eq "SE";
                        $res{RawBases}=$res{RawReads} * $1;
                        $readsLen1=$1;
                }
                if(/^\%GC\s+(\d+)/){
                        $res{GC}=$1;
                }
                @info=split /\s+/;
                if(/#Quality/){
                        while(<I2>){
                                chomp;@info=split /\s+/;
                                last if $info[0]>29;
                                $posNum{20} -= $info[1] if $info[0] <20;
                                $posNum{30} -= $info[1] if $info[0] < 30;
                        }
                        $res{Q20}=sprintf("%.2f",$posNum{20}/$total*100);
                        $res{Q30}=sprintf("%.2f",$posNum{30}/$total*100);
                        last;
                }
        }
        if($ps eq "PE"){
        while(<I3>){
                chomp;
                @info=split /\s+/;

                if(/Total\sSequences\s(\d+)/){
                        warn "paired reads with different reads number\n$fastqcRes\n" if $readsNum1 != $1;
                        $total=$1;
                        $posNum{20}=$total; $posNum{30}=$total;
                        $res{RawReads} += $total;
                }
                if(/Sequence\slength\s(\d+)[-_](\d+)/){
                        warn "can not statistic the base number, because the reads length is not the same in one fastq file\n$fastqcRes\n";
                }
                if(/Sequence\slength\s(\d+)/){
                        warn "paired reads with different reads length\n$fastqcRes\n" if $readsLen1 != $1;
                        $res{RawBases} += $total * $1;
                }

                if(/^\%GC\s+(\d+)/){
                        $res{GC}=$res{GC}."_".$1;
                }
                if(/#Quality/){
                        while(<I3>){
                                chomp;@info=split /\s+/;
                                last if $info[0]>29;
                                $posNum{20} -= $info[1] if $info[0] <20;
                                $posNum{30} -= $info[1] if $info[0] < 30;
                        }
                        $res{Q20}=$res{Q20}."_".sprintf("%.2f",$posNum{20}/$total*100);
                        $res{Q30}=$res{Q30}."_".sprintf("%.2f",$posNum{30}/$total*100);
        #print O "Sample\tInsertSize\tReadLength\tGC%\tRawReads\tRawBases\tQ20%\tQ30%\n";
                        print O "$samp\t$res{InsertSize}\t$res{ReadLength}\t$res{GC}\t$res{RawReads}\t$res{RawBases}\t$res{Q20}\t$res{Q30}\n";
                        last;
                }
        }}
        if($ps eq "SE"){
                print O "$samp\t$res{ReadLength}\t$res{GC}\t$res{RawReads}\t$res{RawBases}\t$res{Q20}\t$res{Q30}\n";
        }

}
close O;
close I;