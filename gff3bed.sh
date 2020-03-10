# download miRBase gff3 (hg38 based) and transfer to bed (hg38 based)

wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
grep -v primary hsa.gff3 | grep -v '#'| awk '{print $1,$4,$5,$9}' OFS="\t" >hsa.gff3.hg38.bed
