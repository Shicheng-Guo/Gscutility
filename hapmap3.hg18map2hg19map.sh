
	wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
  
	plink --bfile hapmap3_r1_b36_fwd_consensus.qc.poly.recode  --recode --tab --out hapmap3_r1_b36_fwd_consensus.qc.poly.recode
	awk '{print "chr"$1,"\t",$4,"\t",$4+1,"\t",$2}' hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map > hapmap3.hg18.bed
	./liftOver hapmap3.hg18.bed hg18ToHg19.over.chain.gz hapmap3.hg19.bed unmap
	awk '{print $1,"\t",$4,"\t",0,"\t",$2}' hapmap3.hg19.bed > hapmap3_r1_b37_fwd_consensus.qc.poly.recode.map
	perl -p -i -e 's/chr//g' hapmap3_r1_b37_fwd_consensus.qc.poly.recode.map
	mv hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped hapmap3_r1_b37_fwd_consensus.qc.poly.recode.ped
	plink --file hapmap3_r1_b37_fwd_consensus.qc.poly.recode --make-bed --out hapmap3_r1_b37_fwd_consensus.qc.poly.recode
	
