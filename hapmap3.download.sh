# https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3
# b36= GRCh36 =hg18. hg19 == GRCh37. 
# each population

cd /home/local/MFLDCLIN/guosa/hpc/db/hapmap3/vcf
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd.qc.poly.tar.bz2
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/relationships_w_pops_051208.txt
tar xjvf hapmap3_r1_b36_fwd.qc.poly.tar.bz2
tar xjvf relationships_w_pops_051208.txt

# consensus
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped.bz2
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2
bzip2 -d hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped.bz2
bzip2 -d hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2

wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/relationships_w_pops_051208.txt
