#!/usr/bin/sh
# Usage: sh bedgraph2amf.sh -b input.bedgraph -d input.bed 
# Extension: perl ~/bin/tab2matrix.pl > matrix.txt
while getopts b:d: option
do
 case "${option}"
 in
 b) bedgraph=${OPTARG};;
 d) bed=${OPTARG};;
 esac
done
# help Alice to treat the methylation bedgraph data to AMF data
sort -k1,1 -k2,2n $bedgraph > $bedgraph.sort
awk '{print $1,$2,$3,$4}' OFS="\t" $bedgraph.sort > $bedgraph.bed4
awk '{print $1,$2,$3,$1":"$2"-"$3}' $bed > $bed.bed4
bedGraphToBigWig $bedgraph.bed4 ~/work/db/mm9/mm9.chrom.sizes $bedgraph.bw
bigWigAverageOverBed $bedgraph.bw $bed.bed4 $bedgraph.tab

