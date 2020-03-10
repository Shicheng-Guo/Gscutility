# download miRBase gff3 and transfer to bed

wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
grep -v primary hsa.gff3 | grep -v '#'| awk '{print $1,$4,$5,$9}' OFS="\t" >hsa.gff3.hg38.bed
