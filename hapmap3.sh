# download database and script
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/GscPythonUtility/master/liftOverPlink.py
# download hapmap3 data in plink format
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped.bz2
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/relationships_w_pops_051208.txt
bzip2 -d hapmap3_r1_b36_fwd_consensus.qc.poly.recode.ped.bz2
bzip2 -d hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2
# convert from hg18 to hg19 plink file
plink --bfile hapmap3_r1_b36_fwd_consensus.qc.poly.recode --recode --out hapmap3.hg18
./liftOverPlink.py -m hapmap3.hg18.map -p hapmap3.hg18.ped -o hapmap3.hg19 -c hg18ToHg19.over.chain.gz -e ./liftOver
./liftOverPlink.py -m hapmap3.hg18.map -p hapmap3.hg18.ped -o hapmap3.hg38 -c hg19ToHg38.over.chain.gz -e ./liftOver

# update plink to binary mode
plink --file hapmap3.hg19 --make-bed --allow-extra-chr --out hapmap3.hg19
plink --file hapmap3.hg38 --make-bed --allow-extra-chr --out hapmap3.hg38

# hapmap3 data cleaning and filtering
plink --bfile hapmap3_r1_b36_fwd_consensus.qc.poly.recode  --missing
plink --file hapmap3_r1_b36_fwd_consensus.qc.poly.recode --maf 0.01 --make-bed --indep 50 5 2 --out hapmap3_r1_b36_fwd_consensus.qc.poly.recode
plink --bfile hapmap3_r1_b36_fwd_consensus.qc.poly.recode --extract hapmap3_r1_b36_fwd_consensus.qc.poly.recode.prune.in --genome --min 0.185
perl ./run-IBD-QC.pl plink
