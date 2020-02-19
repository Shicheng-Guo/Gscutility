for i in `ls gnomad.genomes.r2.1.sites.chr*.rec.hsa.gff3.sort.rmdup.biallelic.vcf.bgz`
do 
bcftools query -f '%CHROM\t%POS\t%ID\n' $i | awk '{print $1,$2-1,$2,$3}' OFS="\t"
done
