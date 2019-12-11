## 12/11/2019
## 660W cytokine data 
## hg18 to hg19
## imputation and phasing
# download database and script
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/GscPythonUtility/master/liftOverPlink.py
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/ibdqc.pl
# rebuild plink file to avoid chromsome-miss-order problem
plink --bfile Th17Set1 --make-bed --out Th17Set1.sort
# space to tab to generate bed files for liftOver from hg18 to hg19
plink --bfile Th17Set1.sort --recode tab --out Th17Set1.tab
# apply liftOverPlink.py to update hg18 to hg19 or hg38
# Only works in BIRC10, Not HPC, Caused by Python version
mkdir liftOver
./liftOverPlink.py -m Th17Set1.tab.map -p  Th17Set1.tab.ped -o Th17Set1.hg19 -c hg18ToHg19.over.chain.gz -e ./liftOver
