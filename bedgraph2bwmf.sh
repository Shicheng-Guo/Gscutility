for i in `ls *.bedGraph`
do 
(awk 'NR!=1' $i| sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4}') > $i.dehsort 
bedGraphToBigWig $i.dehsort ~/oasis/db/hg19.chrom.sizes $i.dehsort.bw 
bigWigAverageOverBed $i.dehsort.bw  /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed  $i..dehsort.bw.mf
rm $i.dehsort
done
