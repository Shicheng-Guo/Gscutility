# build the reference and annotation for the bioinformatics analysis

mkdir $HOME/db

cd $HOME/db
for i in hg18 hg19 hg38 mm9 mm10 
do
mkdir $i
done

# get human and mouse reference

wget 
