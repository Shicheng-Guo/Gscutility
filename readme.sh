cd /home/guosa/hpc/db/hg19/plan2/7mer
GTACGCA.positions.bed
GATCGCA.positions.bed
/home/guosa/hpc/db/hg19/wgbs

cd /home/local/MFLDCLIN/guosa/hpc/db/hg19/plan2
for i in `ls *bed.sort`
do 
awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' $i > $i.bed
done

ssh nu_guos@transfer.chtc.wisc.edu 

for i in chr{1..22} chrX chrY
do
perl dinucleotideFinder.pl ../fa/$i.fa CG
done

for i in chr{1..22} chrX chrY
do
perl ../dinucleotideFinder.pl $i CG  &
done
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.C.positions.txt --counts --out chr22.C &
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.G.positions.txt --counts --out chr22.G &

for i in chr{1..22} chrX chrY
do
perl ../dinucleotideFinder.pl chr22 $i  &
done


for i in chr{1..22} chrX chrY
do
for j in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.$i.*.genotypes.vcf.gz --recode --positions chr$i.$j.vcf.C.positions.txt --counts --out chr22.$j.1 &
done
done


#!/usr/bin/sh

chr="chr1"
for i in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
perl ../dinucleotideFinder.pl $chr $i  &
done

for j in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
wc -l $chr.$j.positions.bed
done


for i in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions $chr.$i.vcf.C.positions.txt --counts --out $chr.$i.1 &
done

for i in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions $chr.$i.vcf.G.positions.txt --counts --out $chr.$i.2 &
done

for j in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
wc -l $chr.$j.1.frq.count
done

for j in AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG
do
wc -l $chr.$j.2.frq.count
done

for chr in chr{1..22} chrX chrY
do
for i in GTACGCA
do
perl 7merMotifFinder.pl $chr $i  &
done
done
cat *GTACGCA.positions.bed >> GTACGCA.positions.txt
bedtools sort -i GTACGCA.positions.txt > GTACGCA.positions.bed

for chr in chr{1..22} chrX chrY
do
for i in GATCGCA
do
perl 7merMotifFinder.pl $chr $i  &
done
done
cat *GATCGCA.positions.bed >> GATCGCA.positions.txt
bedtools sort -i GATCGCA.positions.txt > GATCGCA.positions.bed

/home/guosa/hpc/db/hg19/plan2/7mer







for i in chr{} chrX chrY
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.C.positions.txt --counts --out chr22.C &
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.G.positions.txt --counts --out chr22.G &


/home/guosa/hpc/db/1000Genome/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
vcftools --gzvcf /home/guosa/hpc/db/1000Genome --positions chr22.AT.positions.txt --recode  --counts --out chr22.AT.vcf 

vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --positions test.pos --recode --counts --out chr22

vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.C.positions.txt --counts --out chr22.C &
vcftools --gzvcf /home/guosa/hpc/db/1000Genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --recode --positions chr22.CG.vcf.G.positions.txt --counts --out chr22.G &





for i in `ls *bed.sort.bed.uni.bed`
do 
sort -u $i >> allCpG-SNP.hg19.bed
done

bedtools sort -i allCpG-SNP.hg19.bed > allCpG-SNP.hg19.bed.sort
bedtools cluster -d 1000 -i  allCpG-SNP.hg19.bed.sort >  allCpG-SNP.hg19.G1000.bed.sort.bed
bedtools cluster -d 500 -i  allCpG-SNP.hg19.bed.sort >  allCpG-SNP.hg19.G500.bed.sort.bed
bedtools cluster -d 250 -i  allCpG-SNP.hg19.bed.sort >  allCpG-SNP.hg19.G250.bed.sort.bed &

Rscript --vanilla summary.R allCpG-SNP.hg19.G1000.bed.sort.bed.summary
Rscript --vanilla summary.R allCpG-SNP.hg19.G500.bed.sort.bed.summary
Rscript --vanilla summary.R allCpG-SNP.hg19.G250.bed.sort.bed.summary


* Web server with PHP 7.0.0 or HHVM 3.18.5 or higher.
* A SQL server, the following types are supported
** MySQL 5.5.8 or higher
** PostgreSQL 9.2 or higher
** SQLite 3.3.7 or higher
** Oracle 9.0.1 or higher
** Microsoft SQL Server 2005 (9.00.1399)


http://www.10.103.135.149/index.php

sudo mount -t cifs //luministuscdata.file.core.windows.net/luministphase1data ./luministuscdata -o vers=3.0,username=luministuscdata,password=02grOPN59qcoAQ2EVaSad1z/28YkDe7j1SA6woz36VnbIdXkfhn8tf40JB+WPRvrUnNGWo7SQRTLaJRGANjH1Q==,dir_mode=0777,file_mode=0777,sec=ntlmssp

sudo mount -t cifs //luministuscdata.file.core.windows.net/sequencingdatarun2 ./luministuscdata2 -o vers=3.0,username=luministuscdata,password=02grOPN59qcoAQ2EVaSad1z/28YkDe7j1SA6woz36VnbIdXkfhn8tf40JB+WPRvrUnNGWo7SQRTLaJRGANjH1Q==,dir_mode=0777,file_mode=0777,sec=ntlmssp

sudo mount -t cifs //luministuscdata.file.core.windows.net/sequencingdatarun3 ./luministuscdata3 -o vers=3.0,username=luministuscdata,password=02grOPN59qcoAQ2EVaSad1z/28YkDe7j1SA6woz36VnbIdXkfhn8tf40JB+WPRvrUnNGWo7SQRTLaJRGANjH1Q==,dir_mode=0777,file_mode=0777,sec=ntlmssp


scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/luministuscdata/*.pdf ./


cutadapt -o ./test/Pool_B-tissue_S5_L002_R1_001.fastq.gz -p ./test/Pool_B-tissue_S5_L002_R2_001.fastq.gz  Pool_B-tissue_S5_L002_R1_001.fastq.gz Pool_B-tissue_S5_L002_R2_001.fastq.gz

cutadapt -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC Pool_B-tissue_S5_L002_R1_001.fastq.gz

cutadapt -o ./test/Pool_B-tissue_S5_L002_R1_001.fastq.gz -p ./test/Pool_B-tissue_S5_L002_R2_001.fastq.gz  Pool_B-tissue_S5_L002_R1_001.fastq.gz Pool_B-tissue_S5_L002_R2_001.fastq.gz

perl smartbismark.pl input2.txt submit
qsub Pool_B-tissue_S5_L002_R1_001.fastq.gz.job

 
for i in `ls Pool*job`
do
qsub $i
done

scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/data2/input2.txt ./


-rwxrwxrwx 1 root root  50G Jun  6 20:22 Pool_C-plasma_S3_L001_R2_001.fastq.gz
-rwxrwxrwx 1 root root  53G Jun  6 20:22 Pool_B-plasma_S2_L001_R2_001.fastq.gz
-rwxrwxrwx 1 root root  46G Jun  6 20:22 Pool_B-plasma_S2_L001_R1_001.fastq.gz
-rwxrwxrwx 1 root root  44G Jun  6 20:22 Pool_A-plasma_S1_L001_R2_001.fastq.gz
-rwxrwxrwx 1 root root  44G Jun  6 20:22 Pool_C-plasma_S3_L001_R1_001.fastq.gz
-rwxrwxrwx 1 root root  39G Jun  6 21:29 Pool_A-plasma_S1_L001_R1_001.fastq.gz
-rwxrwxrwx  1 root   root    40G Jun 15 18:44 Pool_A-tissue_S4_L002_R1_001.fastq.gz
-rwxrwxrwx  1 root   root    44G Jun 15 18:44 Pool_A-tissue_S4_L002_R2_001.fastq.gz
-rwxrwxrwx  1 root   root    42G Jun 15 18:44 Pool_B-tissue_S5_L002_R2_001.fastq.gz
-rwxrwxrwx  1 root   root    38G Jun 15 18:44 Pool_B-tissue_S5_L002_R1_001.fastq.gz
-rwxrwxrwx  1 root   root    48G Jun 15 18:44 Pool_C-tissue_S6_L002_R1_001.fastq.gz
-rwxrwxrwx  1 root   root    53G Jun 16 02:51 Pool_C-tissue_S6_L002_R2_001.fastq.gz

hpc/tools/FastQC

ssh 'guosa@mfldclin.org'@10.103.160.225
ssh -L 22:23.99.137.107:22 -p 22 nu_guos:submit-3.chtc.wisc.edu

ssh nu_guos@submit-3.chtc.wisc.edu 'shg047@23.99.137.107 "cp /home/shg047/data1/fastq/input.txt"' > input.txt
scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu 1 nc %h %p' user@remote2:path/to/file .

scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/data1/fastq/* ./
scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/data2/*gz ./
scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/bin/* ./

scp -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:/home/shg047/bin/* ./

fastqc_v0.11.7.zip

 scp nu_guos@submit-3.chtc.wisc.edu:/home/nu_guos/tools/bowtie2-2.3.4.1-linux-x86_64.zip
 scp nu_guos@submit-3.chtc.wisc.edu:/home/nu_guos/tools/*.zip
 scp nu_guos@submit-3.chtc.wisc.edu:/home/nu_guos/tools/*.bz2
 scp nu_guos@submit-3.chtc.wisc.edu:/home/nu_guos/tools/*.gz


scp .ssh/id_rsa.pub -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:~/.ssh/authorized_keys
scp .ssh/id_rsa.pub -o 'ProxyCommand ssh nu_guos@submit-3.chtc.wisc.edu nc %h %p' shg047@23.99.137.107:~/.ssh/authorized_keys2

ssh-keygen -t rsa
ssh nu_guos@submit-3.chtc.wisc.edu mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'cat >> .ssh/authorized_keys'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'cat >> .ssh/authorized_keys2'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'chmod 700 ~/.ssh'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'chmod 640 ~/.ssh/authorized_keys2'
ssh nu_guos@submit-3.chtc.wisc.edu
alias chtc="ssh nu_guos@submit-3.chtc.wisc.edu"

ssh-keygen -t rsa
ssh shg047@23.99.137.107 mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'cat >> .ssh/authorized_keys'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'cat >> .ssh/authorized_keys2'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'chmod 700 ~/.ssh'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'chmod 640 ~/.ssh/authorized_keys2'
ssh shg047@23.99.137.107 
alias azure="shg047@23.99.137.107"

MiSeq, ~25 million reads,  1.5G, 88Gene, 250K, 0.5G, 3 Sample, 65 hours,750$/3=250$/sample
HiSeq 2500, ~400 million reads, 6 G, 88Gene, 250K, 0.5G, 12 Sample, 24 hours,300$/12=25$/sample
HiSeq 4000, ~ 5 billion reads,  1.5T, 250K, 0.5G, 3000 Sample, 84 hours,20000$/3000=12$/sample
HiSeq X Ten, ~ 6 billion reads, 1.8T, 250K, 0.5G, 3600 Sample, 72 hours,12000/3600=3$/sample

https://designstudio-array.illumina.com/#/upload-targets
plink --bfile plink --extract AllCpGSNP150.hg19.RS_SNP.uni.txt --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input

for i in `ls *cpgsnp.bed`
do
bedtools sort -i $i > $i.sort
done

for i in `ls *cpgsnp.bed.sort`
do
awk '{print $1"\t",$}'
done


for i in chr{1..22} chrX chrY
do
cat $i.hg19_cpgsnp.bed.sort >> AllCpGSNP150.hg19.bed
done


rm AllCpGSNP150.hg19.RS_SNP.txt
rm AllCpGSNP150.hg19.RS_SNP.uni.txt
for i in chr{1..22} chrX chrY
do
awk '{print $7}' $i.hg19_cpgsnp.bed.sort >> AllCpGSNP150.hg19.RS_SNP.txt
done
sort -u AllCpGSNP150.hg19.RS_SNP.txt > AllCpGSNP150.hg19.RS_SNP.uni.txt


/home/local/MFLDCLIN/guosa/hpc/pmrp/phase2/RA


a7:b2:da:ad:ca:09:25:e7:80:8a:59:46:90:31:b1:b6 MFLDCLIN\guosa@birc7-lc

ssh-keygen -t rsa
ssh nu_guos@submit-3.chtc.wisc.edu mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'cat >> .ssh/authorized_keys'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'cat >> .ssh/authorized_keys2'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'chmod 700 ~/.ssh'
cat .ssh/id_rsa.pub | ssh nu_guos@submit-3.chtc.wisc.edu 'chmod 640 ~/.ssh/authorized_keys2'
ssh nu_guos@submit-3.chtc.wisc.edu
alias chtc="ssh nu_guos@submit-3.chtc.wisc.edu"

ssh-keygen -t rsa
ssh shg047@23.99.137.107 mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'cat >> .ssh/authorized_keys'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'cat >> .ssh/authorized_keys2'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'chmod 700 ~/.ssh'
cat .ssh/id_rsa.pub | ssh shg047@23.99.137.107  'chmod 640 ~/.ssh/authorized_keys2'
ssh shg047@23.99.137.107 
alias azure="shg047@23.99.137.107"

bgzip	-c	S_Hebbring_Unr.Guo.vcf	>	S_Hebbring_Unr.Guo.vcf.gz
tabix	-p	vcf	S_Hebbring_Unr.Guo.vcf.gz

vcftools	--gzvcf	S_Hebbring_Unr.Guo.vcf.gz	--remove-indels	--recode	--recode-INFO-all	--out	S_Hebbring_Unr.Guo.reindel
bgzip	-c	S_Hebbring_Unr.Guo.reindel.vcf	>	S_Hebbring_Unr.Guo.reindel.vcf.gz
tabix	-p	vcf	S_Hebbring_Unr.Guo.reindel.vcf.gz


plink	--bfile	..	S_Hebbring_Unr.Guo	--keep	parmloss.txt	--chr	23	--allow-no-sex	--recode	--tab	--transpose	--out	6055529-1-0224008846

setwd("C:\\Users\\guosa\\Downloads")
data<-read.table("1176608-1-0238062177.tped",head=F)
het<-apply(data,1,function(x)	sum(!	as.character(x[5])==as.character(x[6])))
plot(het~data$V4,col="red",cex=2,xlab="Chromosome	X",ylab="Heterozygote")



/gpfs/home/guosa/hpc/hemochromatosis/haplotype

/home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype/exomechip_SNV_PASS_BEAGLE_chr6_phased_sel2.map


sel.map


methylation related CpG-SNPs in PMRP Project

for i in `ls *map.map`
do
awk '{print "chr"$2,$5,$7,$6}' $i > $i.newbed
done


HLA-B chr6:31323299-31324,734
MICA chr6:31367561-31383090
NOTCH4 chr6:32163725-32165371


NOTCH4 chr6:31323299-32165371
28493789-33425320

use lib qw(~/hpc/tools/Cwd-Ext-1.06/modulos/share/perl5); # You may need to change this path

/gpfs/home/guosa/hpc/tools/Statistics-R-0.02/modulos/share/perl5

Statistics/R.pm

perl Makefile.PL PREFIX=./modulos
make
make install


for i in chr{1..22} X Y
do
awk '{$1==chr$i}' ../snp150.hg19.txt >>chr$i.vcf.bed
done


for i in chr{1..22} chrX chrY chrM
do
awk -v chr="$i" '$1==chr' ../snp150.ls .txt >> $i.vcf.bed
echo $i
done

for i in chrM
do
awk -v chr="$i" '$1==chr' ../snp150.hg19.txt >> $i.vcf.bed
echo $i
done



'A/G' => 'R',
'C/T' => 'Y',
'A/C' => 'M',
'G/T' => 'K',
'C/G' => 'S',
'A/T' => 'W',
'A/C/T' => 'H',
'C/G/T' => 'B',
'A/C/G' => 'V',
'A/G/T' => 'D',
'A/C/G/T' => 'N',

chr17   80679961        80679962
chr17   80679983        80679984
chr17   80679984        80679985
chr17   80679989        80679990
chr17   80679990        80679991
chr17   80680036        80680037

chr19   93503   93504   rs2353749       -       A       G/T


chr19   1939700 1939701
grep 11138 ../snp150.hg19.txt | grep chr18
grep 93503 ../snp150.hg19.txt | grep chr19
method 1: chr19   668570/4804043
chr19   119331  119332




# http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/index.html
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/HumanOmni2.5-4v1_B-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/HumanOmni25-8v1-2_A1-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/humanomniexpress-24v1-0_a-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/HumanOmniZhongHua-8v1-2_A-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/Immuno_BeadChip_11419691_B-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/InfiniumCoreExome-24v1-1_A-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/InfiniumExome-24v1-1_A1-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/InfiniumPsychArray-24v1-1_A2-b37.Source.strand.zip
wget http://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/Consortium-OncoArray_15047405_A-b37.Source.strand.zip

for i in `ls *.strand`
do
awk '{print "chr"$2,"\t",$3-1,"\t",$3,"\t","chr"$2":"$3}' $i > $i.bed
done

for i in `ls *.strand.bed`
do
bedtools sort -i $i > $i.sort
bedtools intersect -wa -a $i.sort -b ../allSNP150.hg19.cpg-allele-v4.txt > $i.sort.CpGSNP.bed
wc -l $i.sort.CpGSNP.bed
done

90545

for i in `ls *.strand.bed`
do
wc -l $i 
done



for i in `ls *fwd.txt`
do
awk '{print $1}' $i >> Hapmap2
done
sort -u Hapmap2 > Hapmap3
mv 

bedtools intersect -wao -a /home/local/MFLDCLIN/guosa/hpc/db/Hapmap/phase2-3/oncotarget.bed -b hapmapPhaseIIISummary.txt.cpgsnp

bedtools intersect -wao -a /home/local/MFLDCLIN/guosa/hpc/db/Hapmap/phase2-3/oncotarget.bed -b allSNP150.hg19.cpg-allele-v4.txt


rs767677267                 

chr1	114446340	114451307	PTPN22
 
plink --bfile chr1 --recode 12 fastphase --make-bed --chr 1 --extract PTPN22.txt --keep CEU.txt --from-bp 114446340 --to-bp 114451307 --out PTPN22


Potenital methylation Loading(PML)

scp shg047@23.99.137.107:/home/shg047/data2/A09_S* ./


scp shg047@23.99.137.107:/home/shg047/data1/*pl ./

scp shg047@23.99.137.107:/home/shg047/pkg/*zip ./
scp shg047@23.99.137.107:/home/shg047/pkg/*bz2 ./

should_transfer_files = YES
transfer_input_files = us.dat, wi.dat
when_to_transfer_output = ON_EXIT
log = job.log
output = job.out
error = job.err
request_cpus = 1
request_memory = 20MB
request_disk = 20MB




/home/shg047/pkg


set C/G as reference allele

plink --vcf PTPN22.CPG.vcf.gz --reference-allele mylist.txt --recode vcf --real-ref-alleles



https://imputation.sanger.ac.uk/
eagle --vcf PTPN22.vcf --geneticMapFile=/home/local/MFLDCLIN/guosa/hpc/tools/Eagle_v2.4/tables/genetic_map_hg19_withX.txt.gz --outPrefix PTPN22.CPG

https://data.broadinstitute.org/alkesgroup/Eagle/downloads/
	

fastPHASE -F 5000 PTPN22.chr-1.recode.phase.inp
bedtools intersect -wao -a target.bed -b commonSNP_hg19_v1_sort_cluster.bed.summary
bedtools intersect -wao -a target2.bed -b commonSNP_hg19_v1_sort.txt 
 
rs10858022
rs1217397
rs971173
rs28381068
rs3761936
rs114661042
rs1217390
ssh nu_guos@submit-3.chtc.wisc.edu
 
ssh nu_guos@submit-3.chtc.wisc.educd 
/home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype

chr12:4,505,480-4,592,607
plink --file exomechip_SNV_PASS_BEAGLE_chr8_phased_sel2 --make-bed --out exomechip_SNV_PASS_BEAGLE_chr8_phased_sel2

bedtools sort -i commonSNP_hg19_v1.txt > commonSNP_hg19_v1_sort.txt
bedtools cluster -d 1600 -i commonSNP_hg19_v1_sort.txt > commonSNP_hg19_v1_sort_cluster.bed


grep rs12201499 allsnp150.hg19 > test.txt
grep rs6913437 allsnp150.hg19 >> test.txt
grep rs11966502 allsnp150.hg19 >> test.txt
grep rs73716691 allsnp150.hg19 >> test.txt
grep rs757671 allsnp150.hg19 >> test.txt

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/recombRate.txt.gz

awk '{print $2,$3,$4,$5,$7,$8,$10}' hg19.commonsnp150 > hg19_commonsnp150_trim.txt


12269 
gcc --version
libcurl4-openssl-dev
apt-get install libcurl4-openssl-dev
sudo apt-get install texlive
sudo apt-get install texlive-fonts-extra

source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library("GEOquery")

gse <- getGEO("GSE50579", GSEMatrix = TRUE)
show(gse)
filePaths = getGEOSuppFiles("GSE21653")
filePaths
dim(pData(gse[[1]]))
head(pData(gse[[1]])[, 1:3])
df1 <- getGSEDataTables("GSE3494")
lapply(df1, head)

cd /home/local/MFLDCLIN/guosa/hpc/pmrp/merge
# make the first time merge to find out allele need to be filped. 
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --bmerge S_Hebbring_Unr  --out PMRP-Phase1-phase2-Full
# then filp anyone dataset
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip PMRP-Phase1-phase2-Full.missnp --make-bed --out 
# then filp anyone dataset and then merge again. 
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Flip --bmerge S_Hebbring_Unr --out PMRP-Phase1-phase2-Full
# This the remaining non-merged alleles should be indels. run indel2indel.pl to change phase 2 indel mode to phase I. 
perl indel2indel.pl > S_Hebbring_Unr.bim.bim
mv S_Hebbring_Unr.bim.bim S_Hebbring_Unr.bim
# merge again for the last time 
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --bmerge S_Hebbring_Unr  --out PMRP-Phase1-phase2-Full
# Now you get the merge dataset




plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Flip --extract PMRP-Phase1-phase2-Full.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip_INDEL
plink --bfile S_Hebbring_Unr --extract PMRP-Phase1-phase2-Full.missnp --make-bed --out S_Hebbring_Unr_INDEL




   --out PMRP-Phase1-phase2-Full




cp /home/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.bed ./
cp /home/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.bim ./
cp /home/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.fam ./

cp /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/S_Hebbring_Unr.fam ./
cp /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/S_Hebbring_Unr.bim ./
cp /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/S_Hebbring_Unr.bed ./

cp /home/guosa/hpc/pmrp/phase2/S_Hebbring_Unr.Guo.Forward.bed ./
cp /home/guosa/hpc/pmrp/phase2/S_Hebbring_Unr.Guo.Forward.bim ./
cp /home/guosa/hpc/pmrp/phase2/S_Hebbring_Unr.Guo.Forward.fam ./

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield  --list-duplicate-vars --out FinalRelease_QC_20140311_Team1_Marshfield
plink --bfile S_Hebbring_Unr.Guo.Forward  --list-duplicate-vars --out S_Hebbring_Unr.Guo.Forward

perl exm2rs.pl FinalRelease_QC_20140311_Team1_Marshfield.bim 
perl exm2rs.pl  

perl bim2indel.pl FinalRelease_QC_20140311_Team1_Marshfield.bim.bim > FinalRelease_QC_20140311_Team1_Marshfield.bim
perl bim2indel.pl S_Hebbring_Unr.Guo.Forward.bim.bim > S_Hebbring_Unr.Guo.Forward.bim

plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out plink
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip PMRP-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip

plink-merge.missnp

grep Multiple plink.log | awk -F\' '{print $2}' > M.txt


plink --make-bed --vcf FinalRelease_QC_20140311_Team1_Marshfield.vcf.gz --out FinalRelease_QC_20140311_Team1_Marshfield
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --make-bed --remove exclude.txt --out FinalRelease_QC_20140311_Team1_Marshfield_Clean
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Clean --list-duplicate-vars
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Clean --flip PMRP-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Clean.flip

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --make-bed --chr 1 --from-bp 9804693 --to-bp 9819993 
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip PMRP_2018_Phase_1_and_Phase_2.Guo-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out PMRP_2018_Phase_1_and_Phase_2.Guo

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Flip --exclude PMRP_2018_Phase_1_and_Phase_2.Guo-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip_demiss
plink --bfile  S_Hebbring_Unr.Guo.Forward --exclude PMRP_2018_Phase_1_and_Phase_2.Guo-merge.missnp --make-bed --out  S_Hebbring_Unr.Guo.Forward_demiss

plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out plink

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip PMRP_2018_Phase_1_and_Phase_2.Guo-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out PMRP_2018_Phase_1_and_Phase_2.Guo

1161 genomic postion annotation mistakes in Phase 2 dataset and I replace them with the phase I annotation position. 358556 exm id were changed to rs id. 



plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out plink
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out Guo

perl bim2indel.pl FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim > FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim.bim
mv FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim.bim FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out Guo

cp S_Hebbring_Unr.Guo.Forward S_Hebbring_Unr.Guo.Forward.bim
perl bim2indel.pl S_Hebbring_Unr.Guo.Forward.bim > S_Hebbring_Unr.Guo.Forward.bim.bim
mv S_Hebbring_Unr.Guo.Forward.bim.bim S_Hebbring_Unr.Guo.Forward.bim

cp FinalRelease_QC_20140311_Team1_Marshfield_Flip FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim
perl exm2rs.pl FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim 
mv FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim.bim FinalRelease_QC_20140311_Team1_Marshfield_Flip.bim


1, first solve the indel problem in phase 2
2, change exm id to rs id in phase 1
3, unify phase 1 and phase 2 indel genotype
4, 

# phase I data
perl exm2rs.pl FinalRelease_QC_20140311_Team1_Marshfield.bim 
perl bim2indel.pl FinalRelease_QC_20140311_Team1_Marshfield.bim.bim > FinalRelease_QC_20140311_Team1_Marshfield.bim
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out plink
grep Multiple plink.log | awk -F\' '{print $2}' > M.txt
perl ReplaceBimPosition.pl S_Hebbring_Unr.Guo.Forward.bim M.txt FinalRelease_QC_20140311_Team1_Marshfield.bim > S_Hebbring_Unr.Guo.Forward.bim.bim
mv S_Hebbring_Unr.Guo.Forward.bim.bim S_Hebbring_Unr.Guo.Forward.bim
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --flip plink-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Flip




plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield_Flip --make-bed --out Guo




data<-read.table("FinalRelease_QC_20140311_Team1_Marshfield.bim.bim")

for(i in 1:nrow(data)){
if(data[i,2])
}
grep "1KG"   | awk '{ print $2 }' |  xargs grep {} S_Hebbring_Unr.Guo.Forward.bim

awk 'NR==FNR{a[$1];next;}$1 in a' FinalRelease_QC_20140311_Team1_Marshfield.bim S_Hebbring_Unr.Guo.Forward.bim

CoreExome_24v1p2_A1_Anno.Guo.2018.March.short.csv


plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Flip.fam --extract rs7550295.txt --recode --tab --out test
plink --bfile S_Hebbring_Unr.Guo.Forward --extract rs7550295.txt --recode --tab --out test
plink --bfile S_Hebbring_Unr.Guo.Forward --bmerge FinalRelease_QC_20140311_Team1_Marshfield --merge-equal-pos --make-bed --out PMRP_2018_Phase_1_and_Phase_2.Guo

FinalRelease_QC_20140311_Team1_Marshfield.bim

/home/guosa/hpc/pmrp/phase1
check raw proble id
/home/guosa/hpc/pmrp/phase1
/home/guosa/hpc/pmrp/phase1
plink --make-bed --vcf FinalRelease_QC_20140311_Team1_Marshfield.vcf.gz --out FinalRelease_QC_20140311_Team1_Marshfield

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Clean --list-duplicate-vars

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield_Clean --flip PMRP-merge.missnp --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield_Clean.flip

phas1<-read.table("/home/guosa/hpc/pmrp/phase1/myplink.bim")
phas2<-read.table("/home/guosa/hpc/pmrp/phase2/S_Hebbring_Rel.Guo.bim")
head(phas1)
head(phas2)



plink --bfile data1 --bmerge data2.ped data2.map --make-bed --out merge
plink --bfile data1 --bmerge data2.ped data2.map --make-bed --out merge



fam<-read.table("FinalRelease_QC_20140311_Team1_Marshfield_Clean.fam",sep="")
saminfo<-read.table("FinalRelease_QC_Phenotypes_Marshfield_20140224_Team1.txt",head=T,sep="\t")
head(fam)
head(saminfo)
fam[,5]=saminfo[match(fam[,1],saminfo[,1]),9]
write.table(fam,file="FinalRelease_QC_20140311_Team1_Marshfield_Clean.fam2",sep=" ",row.names=F,col.names=F,quote=F)


# phase I data

cd /home/local/MFLDCLIN/guosa/hpc/pmrp/phase1
plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --make-bed --out FinalRelease_QC_20140311_Team1_Marshfield
plink2 --bfile FinalRelease_QC_20140311_Team1_Marshfield --pca approx  --maf 0.05 --memory 40000 --mds-plot --out phase1.pca

--mds-plot

setwd("/home/local/MFLDCLIN/guosa/hpc/pmrp/phase1")
eigenvec<-read.table("phase1.pca.eigenvec",head=F)
saminfo<-read.table("FinalRelease_QC_Phenotypes_Marshfield_20140224_Team1.txt",head=T,sep="\t")
sam<-saminfo[match(as.character(eigenvec[,1]),saminfo$Sample.Name),]

pdf("phase1.pca1-population.pdf")
Legends<-unique(data.frame(Population=sam$PrimaryAncestry_Clean,Col=as.numeric(sam$PrimaryAncestry_Clean)))
plot(eigenvec[,4]~eigenvec[,3],cex=0.55,col=as.numeric(sam$PrimaryAncestry_Clean),pch=as.numeric(sam$PrimaryAncestry_Clean),xlab="PC1",ylab="PC2")
legend("topright",legend=Legends$Population,col=Legends$Col,pch=Legends$Col,cex=0.55)
dev.off()

pdf("phase1.pca2-population.pdf")
Legends<-unique(data.frame(Population=sam$PrimaryAncestry_Clean,Col=as.numeric(sam$PrimaryAncestry_Clean)))
plot(eigenvec[,6]~eigenvec[,5],cex=0.55,col=as.numeric(sam$PrimaryAncestry_Clean),pch=as.numeric(sam$PrimaryAncestry_Clean),xlab="PC3",ylab="PC4")
legend("bottomright",legend=Legends$Population,col=Legends$Col,pch=Legends$Col,cex=0.55)
dev.off()

pdf("phase1.pca9-sex.pdf")
Legends<-unique(data.frame(Population=sam$PrimaryAncestry_Clean,Col=as.numeric(sam$PrimaryAncestry_Clean)))
plot(eigenvec[,12]~eigenvec[,11],cex=0.55,col=as.numeric(sam$PrimaryAncestry_Clean),pch=as.numeric(sam$PrimaryAncestry_Clean),xlab="PC3",ylab="PC4")
legend("bottomright",legend=Legends$Population,col=Legends$Col,pch=Legends$Col,cex=0.55)
dev.off()


data<-read.table("FinalRelease_QC_Phenotypes_Marshfield_20140224_Team1.txt",head=T,sep="\t")
rlt<-apply(data[,23:ncol(data)],2,function(x) summary(glm(data$Sex~x))$coefficients[2,4])
rlt<-apply(data[,23:ncol(data)],2,function(x) summary(glm(data$Age_Clean~x))$coefficients[2,4])
rlt<-apply(data[,23:ncol(data)],2,function(x) summary(glm(as.numeric(data$PrimaryAncestry_Clean)~x))$coefficients[2,4])


plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --make-bed --remove exclude.txt --out FinalRelease_QC_20140311_Team1_Marshfield_Clean




# check PC4< -0.1

Phase I	Eigen-value
PC1	16.6006
PC2	11.9216
PC3	8.30266
PC4	6.48869
PC5	6.23224
PC6	6.07716
PC7	5.91141
PC8	5.66893
PC9	5.40722
PC10 5.22692

HAMP
https://www.ncbi.nlm.nih.gov/nuccore/NC_000019.10?report=genbank&from=35282346&to=35285143
MALSSQIWAACLLLLLLLASLTSGSVFPQQTGQLAELQPQDRAGARASWMPMFQRRRRRDTHFPICIFCCGCCHRSKCGMCCKT
http://useast.ensembl.org/Homo_sapiens/Gene/Compara_Tree?db=core;g=ENSG00000105697;r=19:35280716-35285143;time=1530041560143.143


hg19 genomic region for HLA-B and MICA
chr6:31,232,075-31,391,038

awk '$7==CEU'


cd /home/local/MFLDCLIN/guosa/hpc/db/1000Genome

awk '$7=="CEU"' integrated_call_samples_v2.20130502.ALL.ped > CEU.txt
plink --bfile chr6 --make-bed --keep CEU.txt --maf 0.05 --snps-only --chr 6 --from-bp 31232075 --to-bp 31391038 
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("red", "yellow"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-CEU.pdf")
levelplot(input, main="LD between HLA-B and MICA in CHINA", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()


awk '$7=="CHB"' integrated_call_samples_v2.20130502.ALL.ped > CHINA.txt
awk '$7=="CHS"' integrated_call_samples_v2.20130502.ALL.ped >> CHINA.txt
plink --bfile chr6 --make-bed --keep CHINA.txt --maf 0.05 --snps-only --chr 6 --from-bp 31232075 --to-bp 31391038 
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("red", "yellow"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-CHINA.pdf")
levelplot(input, main="LD between HLA-B and MICA in CHINA", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()

awk '$7=="JPT"' integrated_call_samples_v2.20130502.ALL.ped > JPT.txt
plink --bfile chr6 --make-bed --keep JPT.txt --maf 0.05 --snps-only --chr 6 --from-bp 31232075 --to-bp 31391038 
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out HLAB-MICA.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("red", "yellow"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-JPT.pdf")
levelplot(input, main="LD between HLA-B and MICA in CHINA", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()


chr6:31,980,794-33,074,227



awk '$7=="CHB"' integrated_call_samples_v2.20130502.ALL.ped > CHINA.txt
awk '$7=="CHS"' integrated_call_samples_v2.20130502.ALL.ped >> CHINA.txt
plink --bfile chr6 --make-bed --keep CHINA.txt --maf 0.05 --snps-only --chr 6 --from-bp 31980794 --to-bp 33074227
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out NOTCH-HLADPQ.input
plink --bfile HLAB-MICA.input --r2 --matrix

data<-read.table("plink.ld")
library(lattice)
pal <- colorRampPalette(c("green", "red"), space = "rgb")
input=data.matrix(data)
pdf("LD-matrix-cut4-CHINA-NOTCH.pdf")
levelplot(input, main="LD between NOTCH4 and HLA-D-P-Q in CHINESE", xlab="", ylab="", col.regions=pal(10), cuts=4, at=seq(0,1,0.25))
dev.off()







plink --file fA --merge-list allfiles.txt --make-bed --out mynewdata

10124 Pmrp.Steven.Study.Individual.List.txt


grep rs4475691 /home/local/MFLDCLIN/guosa/hpc/db/hg19/allsnp150.hg19
grep rs4475691 /home/local/MFLDCLIN/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.bim

/home/local/MFLDCLIN/guosa/hpc/pmrp/phase1

my $bim="/home/guosa/hpc/pmrp/phase1/FinalRelease_QC_20140311_Team1_Marshfield.bim";
my $allsnp="/home/guosa/hpc/db/hg19/allsnp150.hg19";

https://watson.hgen.pitt.edu/pub/dbvor/dbvor.1.11.tgz
https://watson.hgen.pitt.edu/pub/slink/apps
https://watson.hgen.pitt.edu/pub/haplo/haplo.HP.Z


6       exm535823       0       32546879        G       A
6       exm535982_ver2  0       32552096        A       G
6       exm2122129_ver4 0       32552121        T       A

for i in `ls *fa`
do
perl cgpositionFinder.pl $i
done


cp /home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype/*newbed

cd /home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype/

plink --file exomechip_SNV_PASS_BEAGLE_chr6_phased_sel2 --make-bed --snps-only --chr 6 --from-bp 32544312 --to-bp 32552716 --out DRB1

awk '{print $2}' plink.bim > mysnps.txt


/home/local/MFLDCLIN/guosa/hpc/db/hg19/fa/chr10.CpG.positions.txt

awk '{print $2,$3,$4,$5}' hg19.commonsnp150 | grep chr10 > chr10.hg19.snp.txt

grep rs35499948 allsnp150.hg19

chr10.hg19.snp.txt



# The easiest way to get lubridate is to install the whole tidyverse:
install.packages("tidyverse")
# Alternatively, install just lubridate:
install.packages("lubridate")
# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("tidyverse/lubridate")


lSLF.AD.Z

# The easiest way to get lubridate is to install the whole tidyverse:
install.packages("tidyverse")
# Alternatively, install just lubridate:
install.packages("lubridate")
# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("tidyverse/lubridate")


install.packages("caret")
install.packages("mlbench")
install.packages("devtools")

library("devtools")
install_github("tidyverse/lubridate")

library(mlbench)
library(caret)
library(ggpubr)
library(car)
library("BBmisc")
shapiro.test(data$GCC.FA.Z)
shapiro.test(powerTransform(data$GCC.FA.Z))
newdata <- preProcess(data[,3:ncol(data)], method=c("BoxCox"))
tidyverse/lubridate

input1<-data$GCC.FA.Z
input2<-powerTransform(data$GCC.FA.Z)
pdf("GCC.FA.Z.pdf")
par(mfrow=c(2,2))
ggqqplot(input1)
ggqqplot(input2)
dev.off()


wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz

phase3.pca.eigenvec

plink --bfile ../S_Hebbring_Unr.Guo --keep parmloss.txt --chr 23 --allow-no-sex --recode --tab --out 1176608-1-0238062177

plink --bfile S_Hebbring_Unr.Guo --snp rs16964862 --allow-no-sex --recode --tab --out rs16964862
plink --bfile S_Hebbring_Unr.Guo.Forward --snp rs16964862 --allow-no-sex --recode --tab --out rs16964862.Forward
plink --bfile S_Hebbring_Unr.Guo --snp exm184366 --allow-no-sex --recode --tab --out exm184366




rs16964862
plink --bfile exomechip_SNV_PASS_BEAGLE_chr13_phased_sel2 --snp rs16964862 --allow-no-sex --recode --tab --out rs16964862
cd /home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype
plink --bfile S_Hebbring_Unr --update-alleles top_to_AB.txt –-make-bed –-out S_Hebbring_Unr.Forward
rs16964862

# hair color
plink --bfile  S_Hebbring_Unr.Guo.Forward --snp rs117322171 --allow-no-sex --recode --tab --out rs117322171
/home/local/MFLDCLIN/guosa/hpc/hemochromatosis/haplotype

grep T,G 

ln -s /mnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data  /home/local/MFLDCLIN/guosa/hpc/pmrp/phase2/rawdata/S_Hebbring_2128_Released_Data

awk '{pirnt $2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7,"\t",$9}' allsnp150

/mqnt/bigdata/Genetic/Projects/S_Hebbring_2128_Released_Data/PLINK_Files/Strand_Translation_Files

# raw plink file of autism 
/mnt/bigdata/Genetic/Projects/Schrodi_Utah_Autism/data

# copy to my folder
cd /home/local/MFLDCLIN/guosa/hpc/autism/data
cp /mnt/bigdata/Genetic/Projects/Schrodi_Utah_Autism/data/All_samples_Exome_QC.bed ./
cp /mnt/bigdata/Genetic/Projects/Schrodi_Utah_Autism/data/All_samples_Exome_QC.fam ./
cp /mnt/bigdata/Genetic/Projects/Schrodi_Utah_Autism/data/All_samples_Exome_QC.bim ./

# change phen xlsx to plink phenotype file
install.packages("openxlsx")

# 314 samples listed in quantitative iffusion tensor MRI data file
# 256 samples have whole-exom sequencing data
# 14 samples don't have MRI quantitative measurement 
# 19 samples were removed since quality control
# 242 samples were included for the assciation study (MRI ~ Allele + age)

# outlier samples
1, 2 low call rate ( u70704cl, u69388s)
2, Gender discrepancy(u38386cl,u65210cl,317814-UW)
3.1, Family(336051,395993 remove, keep 372278)
3.2, Family(370121 remove, keep 386915)
3.3  MZ twin(u28908s remove, keep u28906s-B-Redo)
3.4  duplicated samples (Saliva vs cell line) (u68413d remove, keep u68413s)
4.0  9 PCA outlier(u62997s,u1941001s,u90503s,u59502cl,u64061s,u65457s,u810031s,u810030s,u59504s)
totally, these samples were removed (space sparate): u62997s,u1941001s,u90503s,u59502cl,u64061s,u65457s,u810031s,u810030s,u59504s,u28908s,u68413s,370121,336051,395993,u38386cl,u65210cl,317814-UW,u70704cl,u69388s

Actually, I found the data have already removed 4.0, 3.4, 3.3, 


library("openxlsx")
phen=read.xlsx("DougsPhenotypes_with_avg_age_QTL.xlsx", sheet = 1)
fam<-read.table("All_samples_Exome_QC.fam")
rank1<-match(unlist(lapply(fam[,1],function(x) gsub("u","",strsplit(as.character(x),"-|s|c")[[1]][1]))),phen$ID)
newphen<-cbind(fam,phen[rank1,])
colnames(newphen)[1:2]<-c("FID","IID")
colnames(newphen)[11]<-c("AvgAge")

# prepare covariates plink file (All_samples_Exome_QC.cov, only have age)
cov<-newphen[,c(1,2,11)]
cov[is.na(cov)]<- -9
write.table(cov,file="DougsPhenotypes_with_avg_age_QTL.newphen.cov",sep="\t",quote=F,col.names=c("FID","IID","AvgAge"),row.names=F)

# prepare multiple phenotype plink file (All_samples_Exome_QC.phen)
mphen<-newphen[,c(1,2,13:ncol(newphen))]
mphen[is.na(mphen)]<- -9
write.table(mphen,file="All_samples_Exome_QC.phen",sep="\t",quote=F,col.names=F,row.names=T)

# use perl script to submit job and creat result file with corresponding phenotype names
exc<-read.table("excludeSample.txt")
match(exc[,1],newphen[,1])

# 1031667 variants removed due to MAF<0.01 and 1618874 variants and 249 samples pass filters and QC.  
plink2 --bfile All_samples_Exome_QC --ci 0.95 --genotypic --maf 0.01 --remove excludeSample.txt --allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno 1 --covar All_samples_Exome_QC.cov --linear --out test
plink --bfile All_samples_Exome_QC --ci 0.95 --linear mperm=500  --maf 0.01 --remove excludeSample.txt --allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno 1 --covar All_samples_Exome_QC.cov  --out test
plink --bfile binary_fileset --recode vcf-iid --out new_vcf
plink --bfile All_samples_Exome_QC --ci 0.95 --linear perm --aperm 10 1000000 0.0001 0.01 5 0.001  --maf 0.01 --remove excludeSample.txt --allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno 1 --covar All_samples_Exome_QC.cov  --out test  


plink --bfile All_samples_Exome_QC --snp exm184366 --allow-no-sex --recode --tab --out exm184366
plink --bfile All_samples_Exome_QC --snp kgp8664031 --allow-no-sex --recode --tab --out kgp8664031
plink --bfile All_samples_Exome_QC --snp rs9786510 --allow-no-sex --recode --tab --out rs9786510

plink --bfile All_samples_Exome_QC --snp rs2761764 --allow-no-sex --recode --tab --out rs2761764
plink --bfile All_samples_Exome_QC --snp rs17118281 --allow-no-sex --recode --tab --out rs17118281





awk '$12<0.00005 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt | grep rs6500552 

awk '$12<0.000000031 {print FILENAME,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' *.assoc.linear | grep AD



# find permuation significant assocaiton
0.05/1618874=3.1*10^-8

awk '$4<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm 


awk '$4<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt

awk '$3<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt | grep rs6500552 


kgp8664031

plink2 --bfile All_samples_Exome_QC --impute-sex
 --impute-sex
 
cov<-newphen[,c(1,2,11)]
cov[is.na(cov)]<- -9
mphen<-newphen[,c(1,2,10,13:ncol(newphen))]
mphen[is.na(mphen[,3]),3]<- -9
mphen[mphen[,3]=="Case",3]=1
mphen[mphen[,3]=="Control",3]=0
mphen[is.na(mphen)]<- -9
write.table(mphen,file="All_samples_Exome_QC.phen",sep="\t",quote=F,col.names=T,row.names=F)

# PCA on chrX data
 
7065331-1-0238095238

# collect chrX haplotype and remove missing 
plink --bfile S_Hebbring_Unr.Guo --chr 23 --allow-no-sex --noweb --recode --tab --out chrX
plink --bfile S_Hebbring_Unr.Guo --chr 24 --allow-no-sex --noweb --recode --tab --out chrY

plink --bfile chrXY --pca approx  --maf 0.05 --memory 30000 --out chrXY.pca

plink2 --file chrX --pca approx --maf 0.05 --memory 30000 --out chrX.pca



# IDn plink 23: X, 24: Y, 25: X+Y, 26:M 
# plink --bfile S_Hebbring_Unr.Guo --maf 0.02 --chr 23 --allow-no-sex --noweb --recode --tab --out chrX
# plink --bfile S_Hebbring_Unr.Guo --maf 0.02 --chr 24 --allow-no-sex --noweb --recode --tab --out chrY
system("plink --bfile S_Hebbring_Unr.Guo --maf 0.02 --chr 23 --allow-no-sex --noweb --recode --tab --out chrX")
system("plink --bfile S_Hebbring_Unr.Guo --maf 0.02 --chr 24 --allow-no-sex --noweb --recode --tab --out chrY")
chrx<-read.table("chrX.ped",stringsAsFactors=F,colClasses = c("character"))
chrx<-chrx[,-which(apply(chrx,2,function(x) sum(x==0))>7000)]
chry<-read.table("chrY.ped",stringsAsFactors=F,colClasses = c("character"))
chry<-chry[,-which(apply(chry,2,function(x) sum(x==0))>7000)]
GenderFScore<-read.table("plink.sexcheck",head=T)
saminfo<-read.table("S_Hebbring_Release_Sample_Sheet.txt",head=T,sep="\t")
# chrx het
Tmp<-c()
for(i in 1:nrow(chrx)){
x<-as.character(unlist(chrx[i,seq(7,ncol(chrx),by=2)]))
y<-as.character(unlist(chrx[i,seq(8,ncol(chrx),by=2)]))
Tmp<-c(Tmp,sum(x==y))
print(i)
}
chrXHetRatio=1-Tmp/((ncol(chrx)-6)/2)
Tmp<-unlist(apply(chry[,7:ncol(chry)],1,function(x) sum(x=="0")))
chrYcallrate=1-(Tmp/(ncol(chry)-6))
pedigreeGender=saminfo[match(chrx[,1],saminfo$Sample_Name),]$Gender
Fscore<-GenderFScore[match(chrx[,1],GenderFScore$IID),]$F
rlt<-data.frame(IID=chrx[,1],chrXHetRatio,chrYcallrate,pedigreeGender,Fscore)
write.table(rlt,file="chrXY.gender.prediction.txt",sep="\t",quote=F)

# find permuation significant assocaiton
0.05/1618874=3.1*10^-8

awk '$4<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm 
awk '$4<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt
awk '$3<0.05 {print FILENAME,$1,$2,$3,$4}' *linear.mperm > autism.linear.permuation.sig.txt | grep rs6500552 

plink --file raw-GWA-data --exclude high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --out raw-GWA-data
plink --bfile raw-GWA-data --extract raw-GWA-data.prune.in --genome --out raw-GWA-data
plink --file data --indep 50 5 2

plot(density(Tmp))

saminfo<-read.table("S_Hebbring_Release_Sample_Sheet.txt",head=T,sep="\t")
rlt<-saminfo[match(data[which(Tmp<100),2], saminfo$Sample_Name),]$Gender


mkdir data1
sudo mount -t cifs //discoverydata.file.core.windows.net/discoverydata ./data1 -o vers=3.0,username=discoverydata,password=jyZAwZ5yzqth7jwWSW+XOLCGB4Xsf4NvOwMNF8u7jwymXisTSO6MBT6YligXzmYIC/jsKuUiD84vRsYeoF5IQA==,dir_mode=0777,file_mode=0777,sec=ntlmssp
mkdir data2
sudo mount -t cifs //discoverydata.file.core.windows.net/discoverydata2 ./data2 -o vers=3.0,username=discoverydata,password=jyZAwZ5yzqth7jwWSW+XOLCGB4Xsf4NvOwMNF8u7jwymXisTSO6MBT6YligXzmYIC/jsKuUiD84vRsYeoF5IQA==,dir_mode=0777,file_mode=0777,sec=ntlmssp
mkdir data3
sudo mount -t cifs //discoverydata.file.core.windows.net/discoverydatatissue ./data3 -o vers=3.0,username=discoverydata,password=jyZAwZ5yzqth7jwWSW+XOLCGB4Xsf4NvOwMNF8u7jwymXisTSO6MBT6YligXzmYIC/jsKuUiD84vRsYeoF5IQA==,dir_mode=0777,file_mode=0777,sec=ntlmssp

sudo apt-get install libssl-dev
sudo apt-get install libxml2-dev

install.packages("openssl")
install.packages("git2r")
install.packages("httr")
install.packages("devtools") 
devtools::install_github("Azure/doAzureParallel") 



install pip 
install cutadapt 

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xzvf chromFa.tar.gz
for i in `ls *fa`
do
perl ../../bin/cgpositionFinder.pl $i &
done


for i in `ls *job`
do
sh $i &
done




tail -n 1 *.txt | awk '/^==>/ {a=substr($0, 5, length-8); next} {print a,$1}' | awk '$2>0 {if ($2 !~/chrY/) print $1}' | xargs -I {} rm {}
tail -n 1 *.hap | awk '/^==>/ {a=substr($0, 5, length-8); next} {print a,$1}' | awk '$2>0 {if ($2 !~/chrY/) print $1}' | xargs -I {} rm {}

for i in `ls *hapInfo.txt`
do
#perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/CpGI.hg19.bed4 --output $i.CpGI.hap
#perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/hg19.LINE.bed --output $i.LINE.hap
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/NAS1/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed --output $i.mhb
done
perl ~/bin/hap2mhl.pl > mhl.txt

sudo rstudio-server stop
sudo rstudio-server start
sudo rstudio-server restart

.libPaths( c( .libPaths(), "/media/Home_Raid1/shg047/R/x86_64-pc-linux-gnu-library/3.4") )

# Estellar2016 claim GSM1279519 is a colon tissue, however, in my prediction it is Intestine, therefore, I compare his tissue with roadmap Small_Intestine and Colon
bigWigCorrelate GSM1120321_UCSD.Small_Intestine.Bisulfite-Seq.STL001.wig.gz.bw GSM1010989_UCSD.Sigmoid_Colon.Bisulfite-Seq.STL003.wig.gz.bw                  # 0.9788
bigWigCorrelate GSM1120321_UCSD.Small_Intestine.Bisulfite-Seq.STL001.wig.gz.bw ../../Estellar2016/bw/GSM1279519_CpGcontext.Colon.txt.bedgraph.sort.hg19.bw   # 0.6365
bigWigCorrelate GSM1010989_UCSD.Sigmoid_Colon.Bisulfite-Seq.STL003.wig.gz.bw ../../Estellar2016/bw/GSM1279519_CpGcontext.Colon.txt.bedgraph.sort.hg19.bw     # 0.637193

[[File:Enrichment-of-mhl-ccr.png|200px]]
[[File:Enrichment-mhl-LCP.png|200px]]

# rebuild Rstudio server in genome-miner Version 1.0.136 – © 2009-2016 RStudio, Inc.
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi rstudio-server-1.0.143-amd64.deb

.libPaths( c( .libPaths(), "/media/Home_Raid1/shg047/R/x86_64-pc-linux-gnu-library/3.4") )


Rstudio sever in genome-miner:
http://132.239.25.238:8787/

Rscript ~/bin/mhlpredict.R -f mhl.txt

for i in `ls *hapInfo.txt`
do
#perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/CpGI.hg19.bed4 --output $i.CpGI.hap
#perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/hg19.LINE.bed --output $i.LINE.hap
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/NAS1/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed --output $i.mhb
done
perl ~/bin/hap2mhl.pl > mhl.txt

for i in `ls SRX*hapInfo.txt`
do
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/CpGI.hg19.bed4 --output $i.CpGI.hap
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/work/db/hg19/hg19.LINE.bed --output $i.LINE.hap
perl ~/bin/haptools.pl --input $i --bed /media/Home_Raid1/shg047/NAS1/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed --output $i.mhb.hap
done

grep SINE rmsk.hg19.bed > hg19.SINE.bed
grep LINE rmsk.hg19.bed > hg19.LINE.bed
grep Simple_repeat rmsk.hg19.bed > hg19.LINE.bed

## How to Build your own R package
install.packages("devtools")     # for R package development
install.packages("roxygen2")     # for R package documents (R function => R man)
library("devtools")              # for R package development
library("roxygen2")              # for R package documents (R function => R man)
devtools::document()             # Creat tar.gz source package for distribution
devtools::build()                # Creat tar.gz source package for distribution
install.packages("../monod_1.3.tar.gz")
library("monod")
#' @export                       # before function to export to namespace
data.table::fread()              # Never use library() or require() in a R package!

install.packages("impute")

qsub SRR1035893.job # check what happened!

scp SRP* shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/NAS1/Jenkinson2017NG/samplesheet/
du ./ -h --max-depth 1

cd /media/Home_Raid1/shg047/NAS1/Estellar2016/hapinfo
perl ~/bin/hapinfomergebysrx.pl ../samplesheet/PRJNA229055.txt
for i in `ls SRX*`
do
perl ~/bin/hapinfo2wig.pl $i > $i.bedgraph
echo $i
done
perl ~/bin/

qstat -u shg047 | grep .read | awk '{print $1}' |xargs -I {} qdel {}
qstat -u shg047 | grep trimmed | awk '{print $1}' |xargs -I {} qdel {}
qstat -u shg047 | grep run1 | awk '{print $1}' |xargs -I {} qdel {}
qstat -u shg047 | grep run2 | awk '{print $1}' |xargs -I {} qdel {}

perl ~/bin/hapinfomergebysrx.pl ../samplesheet/SRP072071.txt
perl ~/bin/hapinfomergebysrx.pl ../samplesheet/SRP072075.txt
perl ~/bin/hapinfomergebysrx.pl ../samplesheet/SRP072078.txt
perl ~/bin/hapinfomergebysrx.pl ../samplesheet/SRP072141.txt

wc -l ../samplesheet/SRP072071.txt
wc -l ../samplesheet/SRP072075.txt
wc -l ../samplesheet/SRP072078.txt
wc -l ../samplesheet/SRP072141.txt


grep SRX1651659 ../samplesheet/SRP072071.txt
grep SRX1651659 ../samplesheet/SRP072075.txt
grep SRX1651659 ../samplesheet/SRP072078.txt

scp /media/Home_Raid1/shg047/work/meth450/dyh/20150413/DataResults/1rawdata/idat/3999547166/*idat 


sguo@sph3736:~/Downloads/RA450$ scp shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/work/meth450/dyh/20150514/DataResults/1rawdata/idat/3999423021/*idat ./
+																																																																													
cp /media/Home_Raid1/dinh/RRBS_HaploInfo/Matrices/* /media/Home_Raid1/shg047/NAS1/mhl/

ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =NA,outlier.colour="white")+ coord_flip()


awk 'NR==7258037 {print}' 6-P-1.sorted.clipped.bam.bed

# define RRBS hotspot regions
cd /home/shg047/oasis/monod/bam
for i in `ls *bam`
do
echo $i 
bedtools bamtobed -i $i | awk '{print $1,$2,$3}' OFS="\t" | grep -v 'chrLambdaNEB' | sort -u | sort -k1,1 -k2,2n > $i.bed
bedtools merge -i $i.bed > $i.merge.bed
done


awk '{if ($1>0.1) print $1,$2,$3}' OFS="\t" MM
awk 'i++ {if($1~/ARRS/) print i}' ../../bak/bak.db
awk 'END{if($1 !~/chrY/) print FILENAME,$i,$2,$3,$4}' 
awk 'END{if($1 !~/chrY/) print FILENAME,$i,$2,$3,$4}' *hapInfo.txt
awk -F\| '$3 > 0 { print substr($3,1,6)}' file1
tail -n 1 *.txt | awk '/^==>/ {a=substr($0, 5, length-8); next} {print a,$1}' | awk '$2>0 {if ($2 !~/chrY/) print $1}' | xargs -I {} rm {}


xargs -I {} qsub {}

cov2methcount.sh
#!/bin/bash
# Usage: sh bedgraph2amf.sh -b input.bedgraph -d input.bed 
# Extension: perl ~/bin/tab2matrix.pl > matrix.txt
while getopts i:o: option
do
 case "${option}"
 in
 b) cov=${OPTARG};;
 d) output=${OPTARG};;
 esac
done
awk '{print $1,$2,"+","CpG",$4/100,$5+$6}' $cov >  $output.methcount
pmd -o $output.pmd $output.methcount
hmr -o $output.hmr $output.methcount


cd /media/Home_Raid1/zhl002/NAS1/WGBS/analysis
bedgraph2amf bedgraph bed 
awk -F":|-" '{print ouch$1,$2,$3,$1":"$2"-"$3}' OFS="\t"

#!/bin/bash
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


cd /media/Home_Raid1/zhl002/NAS1/WGBS/analysis
for i in `ls *trim`
do
sort -k1,1 -k2,2n $i > $i.sort
awk '{print $1,$2,$3,$4}' OFS="\t" $i.sort > $i.bed4
bedGraphToBigWig $i.bed4 ~/work/db/mm9/mm9.chrom.sizes $i.bw
bigWigAverageOverBed $i.bw meth_expr_regions.bed.cor.bed $i.tab
perl ~/bin/matrixGeneration.pl > matrix.txt
print $i
done

find ./ -name "*.txt" -mtime +1 -type f | xargs -I {} scp shg047@genome-miner.ucsd.edu:/media/Home_Raid1/zhl002/NAS1/WGBS/analysis unsupervised \;

bedtools 
bedgraphtools
wigtools

methmatrix2gender.R

/home/shg047/oasis/db/FHM.bed


Genomic dataset for Normal Colon tissue
WGBS	SRX332737	SRR949213	Normal Colon
WGBS	SRX332737	SRR949214	Normal Colon
WGBS	SRX332737	SRR949215	Normal Colon
WGBS    SRX381553       GSM1279519      Colon_normal

Genomic dataset for Normal Lung tissue
WGBS      SRX1649893       SRX1649893         Normal Lung    
WGBS      SRX1651655       SRX1651655         Normal Lung
WGBS      SRX1651658       SRX1651658         Normal Lung
WGBS      GSM1279527       SRX381713          Normal Lung

WGBS	  STL001GA-01	   STL001GA-01	      Normal Lung
WGBS	  STL002LG-01	   STL002LG-01	      Normal Lung
WGBS	  N37-Lung	   N37-Lung	      Normal Lung
RRBS	  ENCFF000LVO	   ENCFF000LVO	      Normal Lung
RRBS	  ENCFF000LVR	   ENCFF000LVR	      Normal Lung

SRX1649893
SRX1651655
SRX1651658

awk '{if ($1>0.1) print $1,$2,$3}' OFS="\t" MM

/home/shg047/oasis/Ziller2013/sortbam/hapinfo
SRR949213 SRR949214 SRR949215
cd /media/NAS3_volume2/Dinh/WGBS_LTS33/Hg19/Ziller_Harvard/BAMfiles
scp Colon_Primary_Normal.chr*sorted.clipped.bam* shg047@tscc-login.sdsc.edu:/home/shg047/oasis/Ziller2013/sortbam/SRX332737
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed > saminfo.txt


for i in `ls *bed`; do awk -F"\t|\'|\/" '{print $1,$2,$3,$5/$6,$5,$6}' $i ; done

cd /home/shg047/oasis/Ziller2013/bw
for i in BiSeq_cpgMethylation_BioSam_1121_Colon_Adjacent_Normal.BiSeq.bed GSM1204465_BiSeq_cpgMethylation_BioSam_1120_Colon_Primary_Tumor.BiSeq.bed BiSeq_cpgMethylation_BioSam_157_REMC_19_colonic_mucosa.BiSeq.bed
do
#awk -F"\t|\'|\/" '{print $1,$2,$3,$5/$6}' OFS="\t" $i > $i.bedgraph
#bedGraphToBigWig $i.bedgraph ../../db/hg19/hg19.chrom.sizes $i.bw
bigWigAverageOverBed $i.bw /oasis/tscc/scratch/shg047/monod/hapinfo/Colon77CCP.txt $i.tab
done
perl ~/bin/tab2matrix.pl 



ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1279nnn/GSM1279527/suppl/GSM1279527_CpGcontext.Lung.txt.gz



DATA[match("chr9:103791378-103791447",rownames(DATA)),]
DATA[match("chr12:125051726-125051786",rownames(DATA)),]
DATA[match("chr5:174152050-174152076",rownames(DATA)),]
perl ~/bin/tab2matrix.pl | grep 'chr5:174152050-174152076'


wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1204nnn/GSM1204465/suppl/GSM1204465_BiSeq_cpgMethylation_BioSam_1120_Colon_Primary_Tumor.BiSeq.bed.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1204nnn/GSM1204466/suppl/GSM1204466_BiSeq_cpgMethylation_BioSam_1121_Colon_Adjacent_Normal.BiSeq.bed.gz


cd /oasis/tscc/scratch/shg047/Ziller2013/fastq
perl ~/bin/smartbismark.pl --input saminfo.txt --genome hg19 --server TSCC --submit no --queue hotel
for i in SRR949210 SRR949211 SRR949212 SRR949213 SRR949214 SRR949215
do
qsub $i.pbs
done


/home/shg047/oasis/Chen2016CellResearch/bam
for i in SRR1654403 SRR1654398
do

echo " #!/bin/csh" > $i.job
echo " #PBS -N bam2sort" >> $i.job
echo " #PBS -q hotel" >> $i.job
echo " #PBS -l nodes=1:ppn=8" >> $i.job
echo " #PBS -l walltime=168:00:00" >> $i.job
echo " #PBS -V" >> $i.job
echo " #PBS -M shg047@ucsd.edu" >> $i.job
echo " #PBS -o $i.o" >> $i.job
echo " #PBS -e $i.e" >> $i.job
echo " #PBS -m abe" >> $i.job
echo " #PBS -A k4zhang-group" >> $i.job
echo "cd \$PBS_O_WORKDIR" >> $i.job
echo  samtools sort -@ 8 -o ../sortbam/$i\_bismark_bt2_pe.sort.bam ../bam/$i\_1_val_1_bismark_bt2_pe.nonCG_filtered.bam >> $i.job
echo  samtools index ../sortbam/$i\_bismark_bt2_pe.sort.bam >> $i.job
qsub $i.job
done


SRR1654396_bismark_bt2_pe.sort.bam


# How to add the rsa key of your own linux computer to remove linux server system to avoid input passwd everytime.
# first run the keygen command to generate the keys of your own computer. 
ssh-keygen -t rsa
# tell the remote servers the key of your own computer so that you don't need to input passwd everytime.
cat ~/.ssh/id_rsa.pub | ssh shg047@genome-miner.ucsd.edu 'cat >> ~/.ssh/authorized_keys'
cat ~/.ssh/id_rsa.pub | ssh shg047@tscc-login.sdsc.edu 'cat >> ~/.ssh/authorized_keys'

awk 'i++ {if($1~/ARRS/) print i}' ../../bak/bak.db
find ./ -name ".R" |xargs -I {} grep unsupervised \;
find ./ -name "*.R" | xargs grep 'unsupervised'
find ./ -name "*.R" | xargs grep 'combat'

for i in `ls *hapInfo.txt`
do
perl ~/bin/hapinfo2wig.pl $i > $i.wig
done

perl ~/bin/haptools.pl --input $i --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i 


cd /media/Home_Raid1/shg047/work/monod/hapinfo
perl methHMH.pl 6-T-1.sorted.clipped.bam.hapInfo.txt 6-P-1.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-1 
perl methHMH.pl 6-T-2.sorted.clipped.bam.hapInfo.txt 6-P-2.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-2 
perl methHMH.pl 6-T-3.sorted.clipped.bam.hapInfo.txt 6-P-3.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-3 
perl methHMH.pl 6-T-4.sorted.clipped.bam.hapInfo.txt 6-P-4.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-4 
perl methHMH.pl 6-T-5.sorted.clipped.bam.hapInfo.txt 6-P-5.sorted.clipped.bam.hapInfo.txt excl.list.txt ./caHMH/caHMH-5 
cat caHMH* > HMH.txt
Rscript caHMH.R

#>>>>>>>>
caHMH.R
#>>>>>>>>
data<-read.table("HMH.txt",head=T,sep="\t",as.is=T)
input<-data.matrix(data[,8:ncol(data)])
caHMH<-data[which(apply(input,1,function(x) sum(x)==0)),1]
unique(as.character(caHMH))
write.table(unique(as.character(caHMH)),file="caHMH.rlt.txt",col.names=F,row.names=F,quote=F)

#>>>>>>>>
caHMH.sh
#>>>>>>>>
awk -F":|-" '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t" caHMH.rlt.txt > caHMH.bed
bedtools intersect -wao -a caHMH.bed -b ~/work/db/hg19/hg19_refGene.bed | sort -k1,1n


cor2bed.sh
awk -F":|-" '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t"

chr19	58951755	58951921	chr19:58951755-58951921
scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt ~/work/db/hg19
awk '{print $1,$2,$2+1}' OFS="\t" ~/work/db/hg19/HsGenome19.CpG.positions.txt > HsGenome19.CpG.positions.bed

for i in `ls SRX*hapInfo.txt`
do
if [ -e "$i.hap" ]
then
perl ~/bin/haptools.pl --input $i --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i 
fi
done

cp /home/shg047/oasis/Encode/hapinfo/ENC*hapInfo.txt ../../monod/hapinfo/
cp /home/shg047/oasis/Estellar2016/hapinfo/SRX*.hapInfo.txt ../../monod/hapinfo/

for i in `ls *hapInfo.txt`
do
perl ~/bin/haptools.pl --input $i --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i 
done

for i in `ls N37*hapInfo.txt`
do
perl ~/bin/haptools.pl --input $i --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i 
done

perl ~/bin/haptools.pl --input 6-T-1.sorted.clipped.bam.hapInfo.txt --bed ~/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output test
head 6-T-1.sorted.clipped.bam.hapInfo.txt.hap
head test.hap


/media/12TB_ext/DD_Ext12T/RRBS_MONOD/Plasma_RRBS_2015/BAMfiles

addr:132.249.107.90    2017-05-07

#!/usr/bin/perl
use strict;
chdir "/home/shg047/work/monod/rrbs_kun";
my @file=glob("RRBS*P*");
foreach my $file(@file){
my (undef,$cancerType,$Samid,undef)=split /[-P]/,$file;
my $newName="$cancerType-P-$Samid";
system("cp $file $newName");
}

/media/12TB_ext/DD_Ext12T/Capture_MONOD/150209_SN216/SeqCap/BAMfiles


[59] "WB_centenarian.all_chrs"               
[60] "WB_middle-age.all_chrs"                
[61] "WB_new-born.all_chrs"                  
[62] "SRX381569_tumor_colon"                 
[63] "SRX381716_adenocarcinoma_lung"         
[64] "SRX381719_squamous_cell_tumor_lung"    
[65] "SRX381722_small_cell_tumor_lung" 

find . -name "*" -type d -exec rm -r "{}" \;

for i in `ls *bam.hapInfo.txt`
do
perl ~/bin/haptools.pl --input $i --bed ../mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i
done

wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2600/E-MTAB-2600.sdrf.txt
perl -lane "{print @F[31]}" E-MTAB-2600.sdrf.txt | xargs -I {} wget {}

samtools tview ../sortBam/CTR97_trimmed.fq.gz_bismark_bt2.sort.bam ~/oasis/db/hg19/hg19.fa chr10:100027918-100027944
samtools tview ../bam/6-P-1.sorted.clipped.bam ~/work/db/hg19/hg19.fa chr10:100027918-100027944

cp /media/12TB_ext/DD_Ext12T/RRBS_MONOD/MONOD_RRBS_2014_RemappedTogether_Jan2016/BAMfiles/BAMfiles/*bam ./
cp /media/12TB_ext/DD_Ext12T/RRBS_MONOD/Plasma_RRBS_2016/BAMfiles/*bam ./


perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ../mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed > saminfo.txt

perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt non bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt


for i in `ls *hapInfo.txt`
do
echo "perl methhaptools.pl --input $i --bed ../../mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i"
done


cd /home/shg047/oasis/monod/hapinfo/bak20170501
for i in `ls *hapInfo.txt` 
do
echo " #!/bin/csh" > $i.job
echo " #PBS -N hapinfo2summary" >> $i.job
echo " #PBS -q hotel" >> $i.job
echo " #PBS -l nodes=1:ppn=1" >> $i.job
echo " #PBS -l walltime=72:00:00" >> $i.job
echo " #PBS -V" >> $i.job
echo " #PBS -M shihcheng.guo@gmail.com" >> $i.job
echo " #PBS -m abe" >> $i.job
echo " #PBS -A k4zhang-group" >> $i.job
echo " cd \$PBS_O_WORKDIR" >> $i.job
echo " perl ~/bin/methhaptools.pl --input $i --bed ../../mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed --output $i" >> $i.job
qsub $i.job
done

for i in `ls *hapInfo.txt`
do
echo ""
done

grep -P "chr1" /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.5.bed

awk '{print $1,$2,$3,$1":"$2"-"$3}' /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.5.bed | grep -P "chr1" >  mm9.mhb.0.5.chr1.bed
awk '{print $1,$2,$3,$1":"$2"-"$3}' /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.4.bed | grep -P "chr1" >  mm9.mhb.0.4.chr1.bed
awk '{print $1,$2,$3,$1":"$2"-"$3}' /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.3.bed | grep -P "chr1" >  mm9.mhb.0.3.chr1.bed
perl ~/bin/haptools.pl --input SRR1248477.hapInfo.txt --bed mm9.mhb.0.5.chr1.bed --output SRR1248477.hapInfo.txt.R5
perl ~/bin/haptools.pl --input SRR1248477.hapInfo.txt --bed mm9.mhb.0.4.chr1.bed --output SRR1248477.hapInfo.txt.R4
perl ~/bin/haptools.pl --input SRR1248477.hapInfo.txt --bed mm9.mhb.0.3.chr1.bed --output SRR1248477.hapInfo.txt.R3
bigWigAverageOverBed /home/shg047/oasis/Alice/mouse/hapinfo/Mouse_ESC.meth.bw mm9.mhb.0.5.chr1.bed SRR1248477.hapInfo.txt.R5.tab
bigWigAverageOverBed /home/shg047/oasis/Alice/mouse/hapinfo/Mouse_ESC.meth.bw mm9.mhb.0.3.chr1.bed SRR1248477.hapInfo.txt.R3.tab
bigWigAverageOverBed /home/shg047/oasis/Alice/mouse/hapinfo/Mouse_ESC.meth.bw mm9.mhb.0.4.chr1.bed SRR1248477.hapInfo.txt.R4.tab 


/home/shg047/oasis/Alice/mouse/scMeth/hapinfo

/home/shg047/oasis/Roadmap/tfbs

/media/12TB_ext/DD_Ext12T/Capture_MONOD/141216_HiSeqRapidRun/BAMfiles

for i in `ls *bw`
do 
bigWigAverageOverBed $i ../tfbs/wgEncodeRegTfbsClusteredUnique.bed $i.tab
done


cd /home/shg047/oasis/Roadmap/wig
mkdir mhb
for i in `ls *bw`
do 
bigWigAverageOverBed $i /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed ./mhb/$i.mhb.tab
done

 

R CMD INSTALL openssl_0.9.6.tar.gz --configure-vars='INCLUDE_DIR=/usr/bin/pkg-config LIB_DIR=	'
libssl-dev
system("wget https://cran.r-project.org/src/contrib/xml2_1.1.1.tar.gz")
install.packages("xml2_1.1.1.tar.gz")
system("wget https://cran.r-project.org/src/contrib/selectr_0.3-1.tar.gz")
install.packages("selectr_0.3-1.tar.gz")
system("wget https://cran.r-project.org/src/contrib/magrittr_1.5.tar.gz")
install.packages("magrittr_1.5.tar.gz")
install.packages("rvest_0.3.2.tar.gz")
install.packages("BatchGetSymbols_1.1.tar.gz")

cd /home/shg047/oasis/Jenkinson2017NG/methyfreq


cd /home/shg047/oasis/Jenkinson2017NG/methyfreq
coverage2cytosine --merge_CpG --minDepth 6 --gzip --zero_based --genome_folder /home/shg047/oasis/db/hg19/ -o SRR3263674.test SRR3263674_1_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov.gz


coverage2cytosine --merge_CpG --minDepth --zero_based --gzip --genome_folder $BismarkRefereDb -o SRR3263674.mergeCpG.bed $sample1\_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov.gz



for i in {1..22} X Y M
do
cd /home/shg047/oasis/db/mm9
perl ~/bin/chrosomeCut.pl chr$i 10000 >> mm9.cut10k.bed
done 
wc -l /home/shg047/oasis/db/mm9/mm9.cut10k.bed
cd /home/shg047/oasis/Alice/mouse/scMeth/sortbam
perl /home/shg047/bin/samInfoPrep4Bam2Hapinfo.pl /home/shg047/oasis/db/mm9/mm9.cut10k.bed > saminfo.txt
perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt NON bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/mm9.CpG.positions.txt


# touch all the files 
find . -exec touch {} \;

cp /media/Home_Raid1/shg047/work/Alice/mouse/mhb/mm9.mhb.0.5.bed

for i in `ls *hapInfo.txt`
do
name=${i%%.*}
echo $name $i > $name.input
perl ../hapinfo2mhl.pl $name.input mm9.mhb.0.5.bed > ./MHB/$name.mhb.mhl
done




for i in `ls *hapInfo.txt`
do
name=${i%%.*}
done


scp *.hapInfo.txt shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/NAS1/Alice/mouse/hapinfo/public/
scp SRX209459.hapInfo.txt
scp SRX209458.hapInfo.txt


perl ~/bin/smartMethSRR.pl SRP072071.txt 33 submit
perl ~/bin/smartMethSRR.pl SRP072075.txt 33 submit
perl ~/bin/smartMethSRR.pl SRP072078.txt 33 submit
perl ~/bin/smartMethSRR.pl SRP072141.txt 33 submit

perl ~/bin/smartMethSRR.pl SRP072071.txt 33 none
perl ~/bin/smartMethSRR.pl SRP072075.txt 33 none
perl ~/bin/smartMethSRR.pl SRP072078.txt 33 none
perl ~/bin/smartMethSRR.pl SRP072141.txt 33 none



bedtools shuffle -i /media/Home_Raid1/zhl002/NAS1/WGBS/MHB/Galonska_bivalentdomain.bed -excl

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from  hg19.chromInfo" > hg19.chrom.sizes
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from  mm9.chromInfo" > mm9.chrom.sizes
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from  mm10.chromInfo" > mm10.chrom.sizes

perl hapinfo2mhl.pl input.txt interest2432.bed > output24321.mhl &
perl hapinfo2mhl.pl input.txt Galonska_bivalentdomain.bed > output3004.mhl &
perl methHMHSum.pl input.txt interest2432.bed output24321 &
perl methHMHSum.pl input.txt Galonska_bivalentdomain.bed output3004 &

cd /home/shg047/oasis/Alice/mouse/kun
perl hapinfo2mhl.pl input.txt interest2432.bed > output24321.mhl 

cd /home/shg047/oasis/Alice/mouse/kun
perl hapinfo2mhl.pl input.txt /media/Home_Raid1/zhl002/NAS1/WGBS/MHB/Galonska_bivalentdomain.bed > output3004.mhl &


cron
# sudo vim /etc/crontab
# 17 18 * * * root bash /media/Home_Raid1/shg047/bak/bak/update.sh
# method 1
sudo vim /etc/crontab
37 16   * * *   root    bash /media/Home_Raid1/shg047/bak/bak/bak.sh
# method 2
sudo cp /media/Home_Raid1/shg047/bak/bak/bak.sh /etc/cron.daily/bak
# then restart 
sudo service cron restart
# check cron auto-run list


17 18 * * * root bash /media/Home_Raid1/shg047/bak/bak/update.sh

01 * * * * root echo "This command is run at one min past every hour"
17 8 * * * root echo "This command is run daily at 8:17 am"
17 20 * * * root echo "This command is run daily at 8:17 pm"
00 4 * * 0 root echo "This command is run at 4 am every Sunday"
* 4 * * Sun root echo "So is this"
42 4 1 * * root echo "This command is run 4:42 am every 1st of the month"
01 * 19 07 * root echo "This command is run hourly on the 19th of July"

load("/media/Home_Raid1/zhl002/NAS1/RNA_seq/morehiseq_020617/ipsnt_hiseq_serum.RData")
load("/media/Home_Raid1/zhl002/NAS1/RNA_seq/mips_nt/seurat_allsample_2.RData")
load("/media/Home_Raid1/zhl002/NAS1/RNA_seq/hiseq_020617/seurat_analysis/ipsnt_hiseq_3.RData")
data<-pbmc@scale.data
id1=unlist(lapply(colnames(data),function(x) unlist(strsplit(x,"-"))[2]))
id2=pbmc@ident
group1<-c(which(id2==1),which(id2==2),which(id2==4),which(id2==6))
group2<-c(which(id2==3),which(id2==4),which(id2==7),which(id2==8))
data1=data[,group1]
data2=data[,group2]
Entropy1<-apply(data1,1,function(x) rnashannon(x))
Entropy2<-apply(data2,1,function(x) rnashannon(x))


R2<-read.table("/media/Home_Raid1/zhl002/NAS1/WGBS/ntips.pvalue.R2.txt",head=T,row.names=1)
bed<-data.frame(cor2bed(rownames(R2)),ipsR2[,1],ipsR2[,3])
bed2<-bed2gene(bed,refbed)[,c(1,2,3,4,5,11,12)]
colnames(bed2)<-c("CHR","START","END","iPS-R2","SCNT-R2","Gene","Group")
write.table(bed2,file="LD-R2-table.txt",col.names =NA,row.names =T,sep="\t",quote=F)
LDR2<-read.table("/media/NAS3_volume2/shg047/Alice/mouse/hapinfo/kun/plan3/LD-R2-table.txt")

/media/Home_Raid1/shg047/work/db/mm9/whyte_EN.bed
setwd("/media/NAS1/ZL_NAS1/WGBS")
data=read.table("markers_mhl_regions.txt",head=F)
head(data)
par(mfrow=c(3,1))
boxplot(subset(data,data[,20]=="pluripotent-serum")[,c(4:7)])
boxplot(subset(data,data[,20]=="diff-serum")[,c(4:7)])
boxplot(subset(data,data[,20]=="diff-mix")[,c(4:7)])


source("http://bioconductor.org/biocLite.R")
biocLite("DeconRNASeq")
library("DeconRNASeq")

ES Enhancer

for i in `ls *job`
do
sh $i &
done

perl ~/bin/smartMethSRR.pl SRP072071.txt 33 non
perl ~/bin/smartMethSRR.pl SRP072075.txt 33 non
perl ~/bin/smartMethSRR.pl SRP072078.txt 33 non
perl ~/bin/smartMethSRR.pl SRP072141.txt 33 non

for i in `ls *hapInfo.txt`
do
grep 34303551 $i > ./test/$i.head
done

perl methHMHSum.pl input.txt interest2432.bed output24321
/media/Home_Raid1/zhl002/NAS1/RNA_seq/hiseq_020617/seurat_analysis
/home/shg047/oasis/db/mm9
ipsnt_33k_4.RData
pbmc33k.merged@data
ID=1,2,3,4
indx= 4 5 9 10
SRR3269805_2.fastq.gz

/media/Home_Raid1/shg047/NAS1/Alice/mouse/hapinfo/jpg

perl ~/bin/bam2hapInfo2PBS.pl  saminfo.txt submit bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt

File:CapseqSaminfoConfigBAM2hapinfo2017.txt 

awk '{ sum += $3-$2; n++ } END { if (n > 0) print sum / n; }' CapSeq.bed
samtools view PCT-6.sorted.clipped.bam |awk '{print $3,$4,$4+100}' OFS="\t" | bedtools merge -d 100 -i  - > CapSeqqMerge.bed &

for i in ls *bam
do
samtools view $i | awk "{print $3,$4,$4+100}" 
done

samtools view 


sudo vim /etc/init.d/lampp
#!/bin/bash
/opt/lampp/lampp start
sudo update-rc.d lampp defaults
insserv: warning: script 'K01lampp' missing LSB tags and overrides
insserv: warning: script 'lampp' missing LSB tags and overrides

/media/12TB_ext/DD_Ext12T/Capture_MONOD/141216_HiSeqRapidRun/BAMfiles

cd ..
for i in `ls *bam.hapInfo.txt`
do
head -n 5000 $i > ./test/$i.head
done
cd ./test
perl hapinfo2bed.pl Indx04.sortc.bam.hapInfo.txt.head | sort -k1,1 -k2,2n | bedtools merge -i - > Indx04.bed
bedtools intersect -wa -a ~/work/db/mm9/mm9.refGene.bed -b Indx04.bed > Indx04.input.bed
perl ../methHMHSum.pl input.txt Indx04.input.bed output 


perl methHMHSum.pl input.txt interest12.bed output12-2 &
perl methHMHSum.pl input.txt interest-151.bed output151 &
perl methHMHSum.pl input.txt interest-548.bed output548 &


for i in `ls *hapInfo.txt`
do
grep 34303551 $i > ./test/$i.head
done

34303551

/media/Home_Raid1/shg047/NAS1/Alice/mouse/hapinfo
my ($bin,$NM,$chr,$strand,$txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds,undef,$genesymbol,undef) = split /\t/;
	
chr10:100000000-100010000       CCTCTCT 1       100000149,100000169,100000204,100000238,100000301,100000306,100000361
/media/Home_Raid1/shg047/NAS1/RRBS
scp guoshicheng@gate1.picb.ac.cn:/picb/humpopg6/guoshicheng/lvzhen/*.gz ./

the website:
https://plus.google.com/112774990052239691835
Dongqing Tan
593 6TH AVE
SAN FRANCISCO, CA 94118-3816
United States
Phone number: 415-373-7338

for i in `ls *bam.hapInfo.txt`
do
perl hapinfosummary.pl $i >> result.txt
done

for i in {1..11}
do
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer$i/tracks_hg19/Human_PancreaticCancer$i.meth.bw ./
done

for i in {1..2}
do
http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas$i/tracks_hg19/Human_NormalPancreas$i.meth.bw
done


for i in `ls *bam.hapInfo.txt`
do
head -n 50 $i > ./test/$i.head
done


/media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/public
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDzCsAX+XaSJsDJR71BlYeOLPGgmz5JtpuvL1mYwOsvEZGNA7Y7ufxZE3gDHTZ/fKjfVbewH2FLu+2It6d9LayEy8JQYCcb0oqY1nXHZMdf8bpIfPGz9GJabOHLV6iqx0avpVg64gLtLiWsNV9AfBWd8gjjBTPK33M69SH4/QXAi9zglu6Y325M7FdhyM/7GH/5rIEjd6P5GM6LWJQQkAKpGJf6RcYg9bSyy239SlAFkEJMnBTYCo9xOJYHcyPQAuX6WKiZd/3Sxlqce2zpOIC6O7Q6Q07HWa5aO54R9YXy1pz5NBbn9EtAOZIcS1OuklD2zs8hsjRkijWHlX20Ahbn shg047@tscc-login1.sdsc.edu
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQC33azQGufroerEJp+PuQmrskhbFa6n9sxfXDoKeNI0yPiJCMwvSdxlsZ6k+BqB+84tcTssbHPSVW1EJ5bVRFt6GuDbZoaewHffCkd7hVYZCvZFFtU/o69O+3i0MzWpfO7nwY5Dfhg4rO33Xh8ijKTL4PdAHeoqQ3Cf5aYRfe3f71MEdS7ObfC1FL5h1vKtIsL/IReJ6EUp6h37wyqtnkIQqmkIhJHYc/zCh+CJVgOU6RiCo7rbJBFPozHEITBfrDzhe0azO3+eGdzcerLoKok8eqSc9AKCmSYUW4kdhnuN1j/nOGUaf3i0jZZUWyBldxykxh8Tob1xrgSnL9gA9PfD shg047@genomeMiner

for i in {1..22} X Y
do
grep chr$i HsGenome19.CpG.positions.txt > HsGenome19.CpG.positions.chr$i.txt
done

chr10:116387101-116387201

for i in `ls *hapInfo.txt`
do
perl ~/bin/hapinfo2mf.pl $i > ./methyfreq/$i.bedgraph 
echo $i
done

find . -name "*" -type d -exec rm -r "{}" \;

/media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/public
You can also make a payment via IVR (Interactive Voice Response) by calling 1-800-892-4357 without any fees or charge.

cp /media/Home_Raid1/zhl002/NAS1/projects/mouse_WGBS/bedfiles/mESC_enhancer.bed ~/work/db/mm9
cp /media/Home_Raid1/zhl002/NAS1/projects/mouse_WGBS/bedfiles/Galonska_bivalentdomain.bed ~/work/db/mm9

for i in mESC_enhancer.bed Galonska_bivalentdomain.bed mm9.Enhancer.bed mm9.Promoter.bed mm9.Exon.bed mm9.Intron.bed
do
for j in miPS_B3 miPS_1E12P20 SCNT_NB3.BED SCNT_B12 B6_mESC HA.129.mESC ST_mESC_mix
do
echo $i $j
perl /media/Home_Raid1/shg047/work/Alice/bin/AverageValueNearbyFunctionRegion.pl /media/Home_Raid1/shg047/work/db/mm9/$i  $j.BED.txt.trim  $j.$i.dist.txt
done
done

/media/Home_Raid1/zhl002/NAS1/projects/mouse_WGBS/bedfiles/Galonska_bivalentdomain.bed

miPS_B3.BED.txt.trim

ls miPS_B3.BED.txt.trim
ls miPS_1E12P20.BED.txt.trim
ls SCNT_NB3.BED.txt.trim
ls SCNT_B12.BED.txt.trim
ls B6_mESC.BED.txt.trim
ls HA.129.mESC.BED.txt.trim
ls ST_mESC_mix.BED.txt.trim

echo "# Twin Methylation Dataset Analysis" >> README.md
git init
git add README.md
git commit -m "project start"

/media/Home_Raid1/zhl002/NAS3/ZL_LTS33/mouse_WGBS/bedfiles

for i in `ls *.bed`
do
awk '{if ($4>5) print $1,$2,$3}' OFS="\t" $i > $i.txt
done

#!/usr/bin/perl
use strict;
my @file=glob("*.bed");
foreach my $file(@file){
open F,$file;
while(<F>){
my($chr,$start,$end,$num)=split/\s+/;
print "$chr\t$start\t$end\n" if $num>10;
}
}

/home/shg047/bak/plink/china/gzip/input.bed
cd /home/shg047/bak/plink/china/gzip/
perl -p -i -e 's/\s+/\t/g' /home/shg047/bak/plink/china/gzip/input.bed

# plan 1: 2i vs serum in SCNT
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun
cat Indx06.sortc.bam.hapInfo.txt Indx07.sortc.bam.hapInfo.txt > nt2i.HapInfo.txt
cat Indx09.sortc.bam.hapInfo.txt Indx10.sortc.bam.hapInfo.txt > ntserum.HapInfo.txt
perl ~/bin/hapinfoMerge.pl nt2i.HapInfo.txt > nt2i.HapInfo.Uni.txt
perl ~/bin/hapinfoMerge.pl ntserum.HapInfo.txt > ntserum.HapInfo.Uni.txt
http://www.rencai8.com/job_info?action=view&job_position_id=476884
mv nt2i.HapInfo.Uni.txt ./plan1
mv ntserum.HapInfo.Uni.txt ./plan1
cd ./plan1 
for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl nt2i.HapInfo.Uni.txt > nt2i.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl nt2i.HapInfo.Uni.txt.$i > R2.nt2i.HapInfo.Uni.txt.$i
done

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl ntserum.HapInfo.Uni.txt > ntserum.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ntserum.HapInfo.Uni.txt.$i > ntserum.HapInfo.Uni.txt.$i
done

# plan 2: ips (1,2,4,5) vs scnt (6,7,9,10)
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun
touch  ips.HapInfo.txt
for i in 01 02 04 05 
do
cat Indx$i.sortc.bam.hapInfo.txt >> ips.HapInfo.txt
done

touch  scnt.HapInfo.txt
for i in 06 07 09 10 
do
cat Indx$i.sortc.bam.hapInfo.txt >> scnt.HapInfo.txt
done

perl ~/bin/hapinfoMerge.pl ips.HapInfo.txt > ips.HapInfo.Uni.txt
perl ~/bin/hapinfoMerge.pl scnt.HapInfo.txt > scnt.HapInfo.Uni.txt

mv ips.HapInfo.Uni.txt ./plan2
mv scnt.HapInfo.Uni.txt ./plan2
cd ./plan2

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl ips.HapInfo.Uni.txt > ips.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ips.HapInfo.Uni.txt.$i > R2.ips.HapInfo.Uni.txt.$i
done

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl scnt.HapInfo.Uni.txt > scnt.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl scnt.HapInfo.Uni.txt.$i > R2.scnt.HapInfo.Uni.txt.$i
done

# plan 3: ips(4,5) vs scnt(9,10)
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun
touch ips.HapInfo.txt
for i in 04 05 
do
cat Indx$i.sortc.bam.hapInfo.txt >> ips.HapInfo.txt
done

touch scnt.HapInfo.txt
for i in 09 10 
do
cat Indx$i.sortc.bam.hapInfo.txt >> scnt.HapInfo.txt
done

perl ~/bin/hapinfoMerge.pl ips.HapInfo.txt > ips.HapInfo.Uni.txt
perl ~/bin/hapinfoMerge.pl scnt.HapInfo.txt > scnt.HapInfo.Uni.txt

mv ips.HapInfo.Uni.txt ./plan3
mv ips.HapInfo.Uni.txt ./plan3
cd ./plan3

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl ips.HapInfo.Uni.txt > ips.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ips.HapInfo.Uni.txt.$i > R2.ips.HapInfo.Uni.txt.$i
done

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl scnt.HapInfo.Uni.txt > scnt.HapInfo.Uni.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl scnt.HapInfo.Uni.txt.$i > R2.scnt.HapInfo.Uni.txt.$i
done

cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun/plan1
perl ../../../../bin/R2matrix.pl > plan1.R2.txt
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun/plan2
perl ../../../../bin/R2matrix.pl > plan2.R2.txt
cd /media/Home_Raid1/shg047/work/Alice/mouse/hapinfo/kun/plan3
perl ../../../../bin/R2matrix.pl > plan3.R2.txt


 for i in `ls *.bam`; 
 do 
 touch $i.bam2mf.job
 echo '#!/bin/csh' >$i.bam2mf.job
 echo "#PBS -N $i" >> $i.bam2mf.job
 echo "#PBS -q hotel" >> $i.bam2mf.job
 echo "#PBS -l nodes=1:ppn=1" >> $i.bam2mf.job
 echo "#PBS -l walltime=168:00:00" >> $i.bam2mf.job
 echo "#PBS -M shihcheng.guo@gmail.com" >> $i.bam2mf.job
 echo "#PBS -m abe" >> $i.bam2mf.job
 echo "#PBS -A k4zhang-group" >> $i.bam2mf.job
 echo "cd $(pwd)" >> $i.bam2mf.job
 echo bismark_methylation_extractor --single-end --bedGraph --cutoff 1 --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../sortBam/$i >> $i.bam2mf.job; 
 qsub $i.bam2mf.job
 done
 
perl ~/bin/hapinfo2bedgraph.pl CTR97.hapInfo.txt > CTR97.bedgraph

head /home/shg047/oasis/DennisLo2015/methyfreq/CTR97.merged_CpG_evidence.cov

chr10:100004428-100004429       80.000000
chr10	100004428	100004429	chr10:100004428-100004429

grep 100004429 /home/shg047/oasis/DennisLo2015/methyfreq/CTR97.merged_CpG_evidence.cov
grep 100004428 CTR97.bedgraph

samtools tview ../sortBam/CTR97_trimmed.fq.gz_bismark_bt2.sort.bam ~/oasis/db/hg19/hg19.fa
samtools view -bh ../sortBam/CTR97_trimmed.fq.gz_bismark_bt2.sort.bam chr10:100004428-100004429 -o CTR97.chr10-100004428-100004429.bam
samtools view ../sortBam/CTR97_trimmed.fq.gz_bismark_bt2.sort.bam chr10:100004428-100004429

perl ~/bin/mergedBam2hapInfoV2.pl test.bed CTR97.chr10-100004428-100004429.bam bismark /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt > CTR97.hapInfo.test.txt


/home/shg047/oasis/DennisLo2015/sortBam

/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.5.bed /oasis/tscc/scratch/shg047/Alice/mouse/sortBam/SRX1019866.sortc.bam bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/mm9.CpG.positions.txt > ../hapinfo/SRX1019866.sortc.bam.hapInfo.txt

cd /media/Home_Raid1/shg047/oasis/Alice/mouse/hapinfo/
bedtools intersect -wa -a mm9.mhb.0.5.bed -b /media/Home_Raid1/shg047/work/db/mm9/CpGI.mm9.txt  | wc -l


/media/Home_Raid1/shg047/oasis/Alice/mouse/hapinfo/mm9.mhb.0.5.bed

/home/shg047/oasis/Alice/mouse/mhb/mm9.mhb.0.5.bed


/media/Home_Raid1/shg047/work/db/mm9

表观遗传学科研资讯
crystalqing@yahoo.com
perl hapinfo2mhb.pl -f input.txt -t 0.5 -o mm9.mhb.multinput.bed &
perl ~/bin/hapinfo2mhb.pl mm9.HapInfo.txt.SumUniq 0.5 > mm9.mhb.0.5.bed &
perl ~/bin/hapinfo2mhb.pl mm9.HapInfo.txt.SumUniq 0.3 > mm9.mhb.0.3.bed &


ln -s  /opt/biotools/htseq/bin/htseq-count  htseq-count
ln -s  /opt/biotools/htseq/bin/htseq-qa htseq-qa

/home/shg047/oasis/Alice/mouse/hapinfo


sftp -o Port:8777 'user@domain.com'@domain.com
sftp -o Port:8777 -o User=user@domain.com domain.com
sshpass -p 'Gsc$$8343383' ssh shg047@tscc-login.sdsc.edu

wget --ftp-user='zhang' --ftp-password='HDIIWpP' ftp://169.228.63.66/
wget -r --user='zhang' --password='HDIIWpP' ftp://169.228.63.66/

/home/shg047/software/gmp-4.2.4/build/.libs
/home/shg047/software/gmp-4.2.4
/home/shg047/software/mpfr-2.4.1/build/.libs
/home/shg047/software/mpfr-2.4.1/build

/home/shg047/software/mpfr-3.1.4

export C_INCLUDE_PATH=/home/shg047/software/mpfr-2.4.1/build
export LD_LIBRARY_PATH=/home/shg047/software/mpfr-2.4.1/build/.libs
export LIBRARY_PATH=$LD_LIBRARY_PATH

cd /home/shg047/software/mpc-1.0.3/build
export LD_LIBRARY_PATH=/home/shg047/software/mpfr-2.4.1/build/.libs
../configure --with-gmp-lib=/home/shg047/software/gmp-4.2.4/build/.libs


configure: error: libmpfr not found or uses a different ABI (including static vs shared).


perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ~/db/mm9/mm9.cut10K.bed

cd /home/shg047/oasis/Ziller2013/methyfreq
coverage2cytosine --merge_CpG  --genome_folder  -o 


perl ~/bin/smartbismark.pl --input saminfo.txt --server TSCC --queue hotel --genome mm9
for i in `ls *pbs`; do qsub $i; done
-rw-r--r--  1 shg047 k4zhang-group  1.1G Feb  6 00:30 Indx06.merged.bam_trimmed.fq
-rw-r--r--  1 shg047 k4zhang-group  39G Jan 27 13:13 Indx05.merged.bam.fastq


/home/shg047/oasis/DennisLo2015/methyfreq

coverage2cytosine --merge_CpG --genome_folder ~/oasis/db/hg19/align/bismark/ -o $sample.tmp $file


wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xzvf chromFa.tar.gz
for i in `ls *fa`
do
perl ../../bin/cgpositionFinder.pl $i &
done

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFaMasked.tar.gz

trim_galore --phred33 --fastqc --illumina --rrbs SRR1286404.fastq.gz --output_dir ../fastq_trim
bismark --bowtie2 --phred33-quals --multicore 2 --fastq -L 32 -N 0 -D 5 -R 1 /home/shg047/oasis/db/hg19/align/bismark ../fastq_trim/SRR1286404_trimmed.fq.gz -o ../bam
bismark_methylation_extractor --multicore 3 --single-end --bedGraph --ignore 3 --buffer_size 4G --zero_based --comprehensive --output ../bedgraph  ../bam/SRR1286404_trimmed_bismark_bt2.bam
samtools sort ../bam/SRR1286404_trimmed_bismark_bt2.bam -o ../sortbam/SRR1286404_trimmed_bismark_bt2.sort.bam
samtools index ../sortbam/SRR1286404_trimmed_bismark_bt2.sort.bam

bismark --bowtie2 --phred33-quals --fastq /media/Home_Raid1/shg047/work/db/mm9/bismark  SRX1091397.fastq

cd /media/Home_Raid1/shg047/oasis/Alice/mouse/test

bismark --bowtie2 --phred33-quals --fastq /media/Home_Raid1/shg047/work/db/mm9/bismark  SRX080191.fastq.gz

chr1:98109663-98109763
# check chr1:100006294-100006356 in SRX080191 what happened, why 
# SRX080191.hapInfo.txt:chr1:100000000-100010000  CC      3       100006294,100006294
samtools view SRX080191.sort.bam | less -S 

cd /home/shg047/oasis/Alice/mouse/hapinfo
grep 100006294 SRX080191.hapInfo.txt
samtools view -b -o SRX080191.short.bam SRX080191.sort.bam chr1:100006294-100006356 
samtools tview SRX080191.sort.bam /home/shg047/db/mm9/mm9.fa
samtools tview SRX080191.sort.bam /home/shg047/oasis/db/mm9/mm9.fa
samtools tview Indx01.merged.bam /home/shg047/oasis/db/mm9/mm9.fa


samtools tview SRX080191.sort.bam /home/shg047/db/hg19/hg19.fa
samtools tview SRX080191.sort.bam /home/shg047/oasis/db/hg19/hg19.fa

perl R2matrix.pl > MatrixR2.txt


chr1:100000000-100010000        CC      3       100006294,100006294
chr1:100000000-100010000        CCCC    1       100006294,100006294,100006350,100006350
chr1:100000000-100010000        CCCCCC  1       100006294,100006294,100006350,100006350,100006356,100006356
chr1:100000000-100010000        TT      1       100006294,100006294

samtools tview SRX080191.short.bam ~/work/db/


for i in {1..19} X Y
do
grep -w chr$i Mouse.mm9.HapInfo.txt > Mouse.mm9.chr$i.Hapinfo.txt 
perl ~/bin/hapinfoMerge.pl Mouse.mm9.chr$i.Hapinfo.txt 2> chr$i.err 
perl ~/bin/hapinfo2mhb.pl Mouse.mm9.chr$i.Hapinfo.txt 0.3 > mm9.mhb.chr$i.0.3.txt 
perl ~/bin/hapinfo2mhb.pl Mouse.mm9.chr$i.Hapinfo.txt 0.5 > mm9.mhb.chr$i.0.5.txt  
done

for i in {1..19} X Y
do
grep -w chr$i mm9.HapInfo.txt > Mouse.mm9.chr$i.Hapinfo.txt 
perl ~/bin/hapinfoMerge.pl Mouse.mm9.chr$i.Hapinfo.txt 2> chr$i.err 
perl ~/bin/hapinfo2mhb.pl Mouse.mm9.chr$i.Hapinfo.txt.uni.txt 0.3 > mm9.mhb.chr$i.0.3.txt 
perl ~/bin/hapinfo2mhb.pl Mouse.mm9.chr$i.Hapinfo.txt.uni.txt 0.5 > mm9.mhb.chr$i.0.5.txt  
done


# scRNA-seq could explain methylation haplotype variation

iPS<-read.table("/oasis/tscc/scratch/shg047/Alice/mouse/alice/scRNA/ipsc_R2high.bed")


# 01/27/2017
for i in `ls *pbs`; do qsub $i; done
perl ~/bin/smartbismark.pl --input saminfo.txt --genome mm9 --server TSCC --queue glean
for i in {1..22} X Y M
do
perl ~/bin/genomecut.pl chr$i 10000 /home/shg047/db/mm9/mm9.chrom.sizes > mm9.$i.10k.bed
done
cat mm9.chr*.10k.bed > mm9.chr*.10K.bed

/home/shg047/db/mm9/mm9.cut10K.bed
perl R2matrix.pl > MatrixR2.txt
cat Indx01.hapInfo.txt Indx02.hapInfo.txt Indx04.hapInfo.txt Indx05.hapInfo.txt > iPS.HapInfo.txt
cat Indx06.hapInfo.txt Indx07.hapInfo.txt Indx09.hapInfo.txt Indx10.hapInfo.txt > SCNT.HapInfo.txt
perl ~/bin/hapinfoMerge.pl iPS.HapInfo.txt
perl ~/bin/hapinfoMerge.pl SCNT.HapInfo.txt

mv SCNT.HapInfo.txt.SumUniq SCNT.hapinfo.txt
mv iPS.HapInfo.txt.SumUniq iPS.hapinfo.txt

perl ~/bin/hapinfo2BlocAvgR2.pl iPS.HapInfo.txt.SumUniq > iPS.MHB.R2.txt
perl ~/bin/hapinfo2BlocAvgR2.pl SCNT.HapInfo.txt.SumUniq > SCNT.MHB.R2.txt

# 01/27/2017
for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl SCNT.hapinfo.txt > SCNT.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl SCNT.hapinfo.txt.$i > R2.SCNT.hapinfo.txt.$i
perl ~/bin/randomSampleFromHaploInfo.pl iPS.hapinfo.txt > iPS.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl iPS.hapinfo.txt.$i > R2.iPS.hapinfo.txt.$i
done


# compare two file list in fold A and fold B by diff and ls
diff <(ls ../methyfreq/*cov.gz | awk -F"/" '{print $3}'| awk -F"_" '{print $1}') <(ls *sort.bam | awk -F"_" '{print $1}')

# compare two file list in fold A and fold B by diff and ls
ls ../methyfreq/*cov.gz | awk -F"/" '{print $3}'| awk -F"_" '{print $1}' > a
ls *sort.bam | awk -F"_" '{print $1}' > b
grep -v -f a b 

# compare two file list in fold A and fold B by diff and ls and then qsub job which are not finished. 
grep -v -f a b | awk '{print "qsub "$1"\*.job"}'

diff <(ls ../methyfreq/*cov.gz | awk -F"/" '{print $3}'| awk -F"_" '{print $1}') <(ls *sort.bam | awk -F"_" '{print $1}')
CTR154_trimmed.fq.gz_bismark_bt2.sort.bam

for i in `ls Indx*bam`
do
samtools fastq $i > $i.fastq &
done

gzip $i &

/media/Home_Raid1/dinh/NAS3_volume1_mnt/shg047/GEO/GSE63123_RAW.tar
tar xvf GSE63123_RAW.tar

1, cluster all scRNA by subset genes plu-diff
2, select good cells and separate ips and scnt
3, find differential express genes between ips and scnt
4, show previous knowledge base on this data
5, check Methylation to these DEG
6, find methylation biomarker/matric to show the difference between ips and scnt

/home/shg047/work/DennisLo2015/bam/
 perl ~/software/SNPsplit/SNPsplit --bisulfite --SNP_file CTR151_trimmed.fq.gz_bismark_bt2.sort.bam
samtools view HOT162_trimmed.fq.gz_bismark_bt2.sort.bam | head -n 8000 > ./test/HOT162_trimmed.fq.gz_bismark_bt2.sort.sam
samtools view Pregnancy.8.run1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam | head -n 8000 > ./test/Pregnancy.8.run1.read1_val_1.fq.gz_bismark_bt2_pe.sort.sam

sh CTR84_trimmed.fq.gz_bismark_bt2.sort.bam.bam2mf.job &
sh LTP1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam.bam2mf.job &


CTR132_trimmed.fq.gz_bismark_bt2.sort.sortn.bam

perl ~/bin/bismarkbam2methyfreq.pl --input saminfo.txt --server TSCC --genome hg19 --queue glean 

bismark_methylation_extractor --no_overlap --merge_non_CpG --cutoff 1 --multicore 8 --paired-end --bedGraph --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../bam/HOT215_trimmed.fq.gz_bismark_bt2.sort.bam
perl ~/bin/hapinfo2mhb.pl merge.txt.SumUniq 0.3 > mm9.MHB.0.3.txt   
perl ~/bin/hapinfo2mhb.pl merge.txt.SumUniq 0.5 > mm9.MHB.0.5.txt 
mkdir ../mhb
mv mm9.MHB.0.3.txt ../mhb
mv mm9.MHB.0.5.txt ../mhb/
Added a short File description table in Methylation Haplotype Analysis User Guide.docx
sudo ufw allow 8787
lynx http://<ip>:8787
change chrome to firfox
http:132.239.25.238:8787
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949201/SRR949201_1.fastq.gz 
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949201/SRR949201_2.fastq.gz 
2017-01-17 13:28:06 (0.00 B/s) - "SRR949201_1.fastq.gz" saved [16274079739]
14679170302 Jan 17 14:46 SRR949201_2.fastq.gz
15350187717 Jan 17 15:07 SRR949201_1.fastq.gz
15550948069 Jan 16 18:50 SRR949201_2.fastq.gz
16274079739 Jan 16 18:50 SRR949201_1.fastq.gz
How to prepare bed format for human and mouse refGene

wget -c -O mm9.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz
wget -c -O mm10.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
wget -c -O hg19.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
wget -c -O hg38.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
gunzip -f *.gz
 
fetchChromSizes hg19 > hg19.chom.sizes
fetchChromSizes hg38 > hg38.chrom.sizes
fetchChromSizes mm9 > mm9.chrom.sizes

 
 
#PBS -N Indx04.job
#PBS -q glean
#PBS -l nodes=1:ppn=6
#PBS -l walltime=168:00:00
#PBS -o Indx04.log
#PBS -e Indx04.err
#PBS -V
#PBS -M shicheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group

 
 
 fetchChromSizes mm10 > mm10.chrom.sizes
 
 perl ../mm9/refGene2bed.pl mm9.refGene.txt mm9.chrom.sizes > mm9.refGene.bed
 perl ../mm9/refGene2bed.pl mm10.refGene.txt mm10.chrom.sizes > mm10.refGene.bed
 perl ../mm9/refGene2bed.pl hg19.refGene.txt hg19.chrom.sizes > hg19.refGene.bed
 perl ../mm9/refGene2bed.pl hg38.refGene.txt hg38.chrom.sizes > hg38.refGene.bed
 mv mm9.refGene.bed ../mm9
 cd ../mm9
 
 
 wc -l *.refGene.bed
 
  826320 hg19_refGene.bed
  1372108 hg19.refGene.bed
  1450953 hg38.refGene.bed
   798143 mm10.refGene.bed
   797710 mm9.refGene.bed

  1372104 hg19.refGene.bed
  1449295 hg38.refGene.bed
   798075 mm10.refGene.bed
   797710 mm9.refGene.bed
  4417184 total
 
  1372108 hg19.refGene.bed
  1450953 hg38.refGene.bed
   798143 mm10.refGene.bed
   797710 mm9.refGene.bed
  4418914 total
  
    1372071 hg19.refGene.bed
  1447078 hg38.refGene.bed
   797990 mm10.refGene.bed
   797657 mm9.refGene.bed
  4414796 total

  14679170302 Jan 17 14:46 SRR949201_2.fastq.gz
-rw-r--r-- 1 shg047 k4zhang-group 15350187717 Jan 17 15:07 SRR949201_1.fastq.gz



 wget -c -O mm9.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz
 wget -c -O mm10.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
 wget -c -O hg19.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
 wget -c -O hg38.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

 for i in `ls *.gz`; 
 do 
 SRR=${i%%_*}; 00
 echo $SRR >> list.txt; 
 done
 
 for j in `sort -u list.txt`
 do
 vdb-validate $j
 done
 
 for i in `ls *.gz`; 
 do 
 gunzip -t $i 2 > $i.err &
 done
 find . -name "*gz.err" -type f -size +0c -exec ls -larth {} \;
 

cutadapt: error: In read named 'SRR949201.44244409 D1JR8ACXX130107:7:1108:14746:100581 length=99': length of quality sequence (90) and length of read (99) do not match

less SRR949201_1.fastq.gz |grep -n -A4 'SRR949201.44244409 D1JR8ACXX130107:7:1108:14746:100581 length=99' 
	
perl ~/bin/smartbismark.pl --input PRJNA201480.txt --submit no --genome hg19 --server TSCC
qsub SRR949201.pbs
qsub SRR949202.pbs
qsub SRR949210.pbs
qsub SRR949211.pbs

sh SRR949201.job &
sh SRR949202.job &
sh SRR949210.job &
sh SRR949211.job &

for i in `ls *.gz`
do
md5sum $i >> md5sum.txt
md5sum $i >> md5sum.txt
done

SRR949202_1 missed in fastq_trim folder
_2 missed in fastq_trim folder
SRR949210_1 missed in fastq_trim folder
_2 missed in fastq_trim folder
SRR949211_1 missed in fastq_trim folder
SRR949211_2 missed in fastq_trim folder


coverage2cytosine --merge_CpG SRR949207_1_val_1_bismark_bt2_pe.nonCG_filtered.bismark.cov.gz --genome_folder ~/oasis/db/hg19/align/bismark/ -o SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph

cd /home/shg047/software/R-3.3.2
./configure --prefix=/home/shg047/software/R-3.3.2 '--with-cairo' \
 '--with-jpeglib' '--with-readline' '--with-tcltk' '--with-x=no'\
 '--with-blas' '--with-lapack' '--enable-R-profiling' '--with-tiff=yes'\
 '--enable-R-shlib'\
 '--enable-memory-profiling'
 make clean
 make
 
 
Name is names what the name 
 
Clear idea, clean story, base on truth and get some another truth. (drugs and interaction)
Our value is fix some gap in the academic fileds.
We would better do some thing beyond context dependent. how to interpret the result or study with context dependent. 
 => human mutation. 
1, small group and train and fight yourself
2, train yourself to find the bar.
3, you should know how to fight alone.  
4, build your own contribution. also you need your project huge. that means you need more energy to have some light on
5, people or paper or boht. 

cd /media/Home_Raid1/shg047/software/InfiniumPurify
for i in `ls jhu*`
do
python InfiniumPurify.py -f $i -c LUAD 
done


Bowtie with parameters "-q --phred33-quals -n 1 -e 99999999 -l 25 -I 1 -X 2000 -a -m 15 -S -p 6",
E(i,j)=log2(TPMi,j/10+1), where TPMi,j refers to transcript-per-million (TPM) for gene i in sample j, as calculated by RSEM


head -n 20 SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov
head -n 20 SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph.merged_CpG_evidence.cov


wc -l  SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph.merged_CpG_evidence.cov


1731880
15609274

grep chrUn_gl000214 SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov
grep SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph.merged_CpG_evidence.cov

head SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov
head SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph.merged_CpG_evidence.cov

system("cat $cov >> $SRX.bedgraph");
/home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt

 #!/bin/csh
 #PBS -N MergeCOVbySRX
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group 
cd /home/shg047/oasis/Estellar2016/methyfreq
coverage2cytosine --merge_CpG SRR1035731_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.head --genome_folder ~/oasis/db/hg19/align/bismark/ -o SRR1035809_1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.bedgraph


Chase Sapphire Preferred (CSP)
Amex Centurion
Citi Prestige
Amex Premier Rewards Gold (PRG)
Amex Platinum
Chase Ritz Carlton
Chase Palladium
Amex SPG
Discover it


 perl ~/bin/samInfoPrep4Bam2Hapinfo.pl /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | grep -v PPP3CA
 
qsub SRR949205.pbs
qsub SRR949212.pbs
qsub SRR949215.pbs

curl "https://gdc-api.nci.nih.gov/legacy/files/4a8ffe0d-d7e6-4712-ad04-472955c84c77?fields=cases.samples.portions.analytes.aliquots.submitter_id,cases.samples.sample_type&format=tsv"
curl "https://gdc-api.nci.nih.gov/legacy/files/087ec4fb-a621-4fcf-8276-1c74782bcc2c?fields=cases.samples.portions.analytes.aliquots.submitter_id,cases.samples.sample_type&format=tsv"
 
du ./ -h --max-depth 1

# 2016-12-30
cd /media/Home_Raid1/shg047/work/Roadmap/bw
for i in GSM1010981_UCSD.Adrenal_Gland.Bisulfite-Seq.STL003.wig.gz.bw GSM983649_UCSD.Esophagus.Bisulfite-Seq.STL003.wig.gz.bw GSM1120336_UCSD.Right_Ventricle.Bisulfite-Seq.STL003.wig.gz.bw GSM1120337_UCSD.Right_Ventricle.Bisulfite-Seq.STL003.wig.gz.bw
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/CpGSF.hg19.sort.bed4 $i.tab
done

cd /media/Home_Raid1/shg047/work/Chen2016CellResearch/bw
for i in `ls *bw`
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/CpGSF.hg19.sort.bed4 $i.tab
mv $i.tab /media/Home_Raid1/shg047/work/Roadmap/bw
done

cd /media/Home_Raid1/shg047/work/Roadmap/bw
perl ~/bin/bigWigAverageOverBed2Matrix.pl > CpGI.txt

data<-read.table("CpGI.txt")
group<-unlist(lapply(rownames(data),function(x) unlist(strsplit(as.character(x),":"))[3]))
table(group)
input<-data.frame(Pos=rownames(data),data,group)
head(input)
library("reshape2")
colnames(input)<-c("POS","Kidney","Right_Ventricle 1","Right_Ventricle 2","Kidney_Tumor","Kidney_Normal","Kidney_Tumor","Kidney_Normal","Esophagus","group")
input.long<-melt(input, id.vars=c("group","POS"))
library(ggplot2)
png("rmsk.CpGI.methylation.png")
ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =NA)+ coord_flip()
dev.off()
tapply(input.long$value,input.long$group,function(x) median(x,na.rm=T))
tapply(input.long$value,input.long$variable,function(x) mean(x,na.rm=T))
head(subset(input.long,group=="CpGI"))



# 2016-12-30
cd /media/Home_Raid1/shg047/work/Roadmap/bw
for i in GSM1010981_UCSD.Adrenal_Gland.Bisulfite-Seq.STL003.wig.gz.bw GSM983649_UCSD.Esophagus.Bisulfite-Seq.STL003.wig.gz.bw GSM1120336_UCSD.Right_Ventricle.Bisulfite-Seq.STL003.wig.gz.bw GSM1120337_UCSD.Right_Ventricle.Bisulfite-Seq.STL003.wig.gz.bw
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed $i.tab
done

cd /media/Home_Raid1/shg047/work/Chen2016CellResearch/bw
for i in `ls *bw`
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed $i.tab
mv $i.tab /media/Home_Raid1/shg047/work/Roadmap/bw
done

cd /media/Home_Raid1/shg047/work/Roadmap/bw
perl ~/bin/bigWigAverageOverBed2Matrix.pl > Repeat.txt

data<-read.table("repeat.txt")
group<-unlist(lapply(rownames(data),function(x) unlist(strsplit(as.character(x),":"))[3]))
table(group)
input<-data.frame(Pos=rownames(data),data,group)
head(input)
library("reshape2")
colnames(input)<-c("POS","Kidney","Right_Ventricle 1","Right_Ventricle 2","Kidney_Tumor","Kidney_Normal","Kidney_Tumor","Kidney_Normal","Esophagus","group")
input.long<-melt(input, id.vars=c("group","POS"))
library(ggplot2)
png("rmsk.repeat.methylation.png")
ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =NA)+ coord_flip()
dev.off()
tapply(input.long$value,input.long$group,function(x) median(x,na.rm=T))
tapply(input.long$value,input.long$variable,function(x) mean(x,na.rm=T))
head(subset(input.long,group=="CpGI"))


# 2016-12-30
cd /media/Home_Raid1/shg047/work/db/hg19
perl rmsk.pl -i rmsk.hg19 -o rmsk.hg19.bed
sort -u rmsk.hg19.bed > rmsk.hg19.bed.temp
sort -k1,1 -k2,2n rmsk.hg19.bed.temp > rmsk.hg19.bed
rm rmsk.hg19.bed.temp
awk '{print $12}' rmsk.hg19 | sort -u

cd /media/Home_Raid1/shg047/work/Chen2016CellResearch/bw
for i in `ls *bw`
do
bigWigAverageOverBed $i /media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed $i.tab
done

perl ~/bin/ Repeat.txt

data<-read.table("Repeat.txt")
info<-read.table("/media/Home_Raid1/shg047/work/db/hg19/rmsk.hg19.bed")
group<-unlist(lapply(info[,4],function(x) unlist(strsplit(as.character(x),":"))[3]))
table(group)
input<-data.frame(data,group)

library("reshape2")
colnames(input)<-c("T","N","T","N","group")
input.long<-melt(input, id.vars=c("group"))

library(ggplot2)
pdf("rmsk.methylation.pdf")
ggplot(aes(y = value, x = group, fill = variable, dodge=variable), data = input.long) + geom_boxplot(outlier.shape =NA,outlier.colour="white")+ coord_flip()
dev.off()
##

git pull
git status
git rebase --continue

#

http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://132.239.25.238/shg047/myhub.txt
 
ln -s /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/bw/GSM1546663_5mC-P1-T-corrected.txt.bw

/media/Home_Raid1/shg047/work/h
UCSC Over
mv /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/bw
myhub.txt

cd /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/
for i in `ls *bw`
do
ln -s /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/$i /media/Home_Raid1/shg047/work/hub/hg19/$i
done

sudo ln -s /media/Home_Raid1/shg047/work/Chen2016CellResearch/GSE63183/GSM1546663_5mC-P1-T-corrected.txt.bw /media/Home_Raid1/shg047/work/hub/hg19/GSM1546663_5mC-P1-T-corrected.txt.bw



# 2016-12-29
 liftOver GSE17972.hg18.bedgraph ~/work/db/hg18/hg18ToHg19.over.chain GSE17972.hg19.bedgraph tmp
 sort -u -k1,1 -k2,2n GSE17972.hg19.bedgraph > GSE17972.hg19.bedgraph.sort
 bedGraphToBigWig GSE17972.hg19.bedgraph.sort ~/work/db/hg19/hg19.chrom.sizes GSE17972.hg19.bw
 track type=bigWig color=0,0,255 visibility=2 maxHeightPixels=128:30:11 smoothingWindow=16 windowingFunction=mean name="PBMC" description="Yanhuang-methylome" bigDataUrl=http://132.239.25.238/shg047/NAS3/shg047/Yanhuang2010/GSE17972.hg19.bw

 liftOver GSE17972.hg18.bedgraph ~/work/db/hg18/hg18ToHg38.over.chain GSE17972.hg38.bedgraph tmp
 sort -u -k1,1 -k2,2n GSE17972.hg38.bedgraph > GSE17972.hg38.bedgraph.sort
 bedGraphToBigWig GSE17972.hg38.bedgraph.sort /media/Home_Raid1/shg047/work/db/hg38/hg38.chrom.sizes GSE17972.hg38.bw
 track type=bigWig color=0,0,255 visibility=2 maxHeightPixels=128:30:11 smoothingWindow=16 windowingFunction=mean name="Yanhuang-methylome" description="PBMC" bigDataUrl=http://132.239.25.238/shg047/NAS3/shg047/Yanhuang2010/GSE17972.hg38.bw
 
 GSE17972.hg19.bw
 GSE17972.hg19.bedgraph.sort
 
for i in `ls *job`
do

if [! -e ]

for i in SRR1035834 SRR1035835 SRR1035844 SRR1035831 SRR1035784 SRR1035893 SRR1035884 SRR1035843 SRR1035895 SRR1035832 SRR1035882 SRR1035881 SRR1035845 SRR1035833
do
qsub $i.fastq.download.job
done

for i in `ls *bam`
do
touch $i
done


 install CPAN
 reload cpan

 # qmap to bedgraph, liftOver hg18 to hg19 and hg38, finally sort the bedgraph 

 for i in `ls *.txt`
 do
 echo $i
 perl qmap2bedgraph.pl $i > $i.bedgraph
 done
 
 # merge liftOver 
 for i in `ls *.bedgraph`
 do
 echo $i
 liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.hg19 tmp
 liftOver $i /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain $i.hg38 tmp
 done
 
 # Sort liftOver bedgraph
 for i in `ls *.bedgraph.hg19`
 do
 echo $i
 sort -k1,1 -k2,2n $i > $i.sort
 done

 for i in `ls *.bedgraph.hg38`
 do
 echo $i
 sort -k1,1 -k2,2n $i > $i.sort
 done

 # merge bedgraph by chrosome
 cat *hg18.sort > GSE17972.YanHuang.hg18.bedgraph 
 cat *hg19.sort > GSE17972.YanHuang.hg19.bedgraph 
 cat *hg38.sort > GSE17972.YanHuang.hg38.bedgraph 

 # bedgraph to bigwig 
 for i in `ls *.bedgraph`
 do
 echo $i
 bedGraphToBigWig $i ~/work/db/hg18/hg18.chrom.sizes $i.bw
 bedGraphToBigWig $i ~/work/db/hg19/hg19.chrom.sizes $i.bw
 bedGraphToBigWig $i ~/work/db/hg38/hg38.chrom.sizes $i.bw
 done

 
 
 
 
 for i in `ls *chr10*txt`
 do
 sort -k1,1 -k2,2n $i.bedgraph.hg19 > $i.bedgraph.hg19.sort
 bedGraphToBigWig $i.bedgraph.hg19 ~/work/db/hg19/hg19.chrom.sizes $i.hg19.bw
 done



 
liftOver GSE17972_HUMtg5lib.qmap.chr10.txt.bedgraph /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain GSE17972_HUMtg5lib.qmap.chr10.txt.bedgraph.hg19 tmp
liftOver GSE17972_HUMtg5lib.qmap.chr10.txt.bedgraph /media/Home_Raid1/shg047/work/db/hg18/hg18ToHg19.over.chain GSE17972_HUMtg5lib.qmap.chr10.txt.bedgraph.hg38 tmp


http://132.239.25.238/shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546666_5mC-P2-N-corrected.txt.bw

shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546663_5mC-P1-T-corrected.txt.bw
shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546664_5mC-P1-N-corrected.txt.bw
shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546665_5mC-P2-T-corrected.txt.bw
shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546666_5mC-P2-N-corrected.txt.bw

http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl= 
http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubUrl=http://132.239.25.238/shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/myhub.txt

track type=bigWig name="XX" description="YY" bigDataUrl=http://132.239.25.238/shg047/NAS3/shg047/Chen2016CellResearch/GSE63183/bw/hg19/GSM1546663_5mC-P1-T-corrected.txt.bw


for i in `ls *5mC-P*txt`
do
sort -k1,1 -k2,2n $i >$i.sort
perl tobedgraph.pl $i.sort 
bedGraphToBigWig $i.sort.bedgraph ~/db/hg19/hg19.chrom.sizes $i.bw 
rm $i.sort 
rm $i.sort.bedgraph
done


for i in `ls *.job`
do
sh $i &
done



GetOptions ( "input=s"   => \$input,           # string
             "submit=s"  => \$submit,          # flag
             "genome=s" => \$genome,          # string
             "server=s" => \$server)          # flag
o
perl ~/bin/smartbismark.pl --input PRJNA201480.txt --submit no --genome hg19 --server TSCC

bismark --bowtie2 --multicore 4 --phred33-quals --fastq -L 25 -N 1 /home/shg047/db/hg19/bismark/ -1 ../fastq_trim/SRR949204_1_val_1.fq.gz -2 ../fastq_trim/SRR949204_2_val_2.fq.gz -o ../bam &


Memo for today's discussion: 
Human brain includes two broad classes of cells: neurons and glial cells. 
New understanding to AS based on deconvolution:
Region projection for AD based on region reference => which regions => Glial cells  or Neurons cells? 8224138
Glial cells  => which subtype ? => inflammatory related pathway? => eQTL => AD
Neurons cells =>  which subtype ? => which pathways ? => eQTL => AD
New Method based on deconvolution:
1,  Brain region (reference) and cell type (reference) projection for AD
2, Differential gene expression after deconvolution identify 
3, Causal network after deconvolution identify cell type (reference) specific interaction? 
Thanks. 

Shicheng


cp /media/Home_Raid1/zhl002/NAS3/WGBS/permutation/*bed  /opt/lampp/htdocs/


cd /opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/hg19
for i in `ls *bed`
do
genome="hg19"
wget -O $i.results.tsv "http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=$genome&requestName=Example+Data&requestSender=Client+A&requestURL=http%3A%2F%2F132.239.25.238%2Fshg047%2FNAS3%2FAlice%2FWGBS%2Fpermutation%2Fhg19%2F$i"
done

cd /opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/mm9
for i in `ls *bed`
do
genome="mm9"
wget -O $i.results.tsv "http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=$genome&requestName=Example+Data&requestSender=Client+A&requestURL=http%3A%2F%2F132.239.25.238%2Fshg047%2FNAS3%2FAlice%2FWGBS%2Fpermutation%2Fmm9%2F$i"
done


/opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/hg19

/opt/lampp/htdocs/shg047/NAS3/Alice/WGBS/permutation/

报案：受到台湾公民（蘇勝慧）的骚扰，诽谤，损害名誉权，人身健康威胁及恐吓案件
chr14:29318659
fastq-dump --maxSpotId 100 --minSpotId 200 SRR1648428

samtools tview N22_bismark_bt2_pe.sort.bam ~/db/hg19/

qstat -u shg047 | grep condo | awk '{print $1}' | xargs -I {} qdel {}﻿
for i in 0 T U  
do 
rm *$i.bam
done

trim_galore --paired --phred33 --fastqc --illumina N21_R1.fastq.gz N21_R2.fastq.gz --output_dir ../fastq_trim
bismark --bowtie2 --multicore 4 --phred33-quals --fastq -L 21 -N 1 ~/db/hg19/align/bismark/ -1 ../fastq_trim/N21_R1_val_1.fq.gz -2 ../fastq_trim/N21_R2_val_2.fq.gz -o ../bam
filter_non_conversion --paired ../bam/N21_R1_val_1_bismark_bt2_pe.bam
samtools sort ../bam/N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bam ../sortbam/N21_bismark_bt2_pe.sort
samtools index ../sortbam/N21_bismark_bt2_pe.sort.bam
bismark_methylation_extractor --cutoff 10 --paired-end --bedGraph --ignore 1 --buffer_size 4G --zero_based --comprehensive --output ../methyfreq  ../bam/N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bam
bedtools intersect -wa -a N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov -b ../bed/target.bed > ../bedgraph/N21.bedgraph

N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov
N21_R1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov.

chr8:67344678-67344679

cd /home/shg047/work/Minghua2016/methyfreq
for i in `ls *bismark.zero.cov`
do
bedtools intersect -wa -a $i -b ../target.bed > $i.bedgraph
done

bedtools intersect -wa -a T9_R1_val_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.cov -b ../target.bed

trim_galore --paired --phred33 --fastqc --illumina SRR949193_1.fastq.gz SRR949193_2.fastq.gz --output_dir ../../fastq_tri
CU1_R1_val_1_bismark_bt2_pe.nonCG_filtered.bam
CU1_val_1_bismark_bt2_pe.nonCG_filtered.bam
trim_galore --paired --phred33 --fastqc --illumina SRR949197_1.fastq.gz SRR949197_2.fastq.gz --output_dir ../fastq_trim
SCNT.wnt5a_mouse.sort.bam
SRR949205_1_val_1_bismark_bt2_pe.bam

/media/Home_Raid1/shg047/db/aligndb/hg19/bismark
scp SRR949206_1_trimmed.fq.gz shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/work/Ziller2013/fastq_trim 
scp SRR949206_2_trimmed.fq.gz shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/work/Ziller2013/fastq_trim
scp SRR949206.bismark.job shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/work/Ziller2013/fastq


zcat ../fastq_trim/SRR949205_1_val_1.fq.gz | head -n 4000 > ../SRR949205_1_val_1.fq
zcat ../fastq_trim/SRR949205_2_val_2.fq.gz | head -n 4000 > ../SRR949205_2_val_2.fq

gzip ../SRR949205_1_val_1.fq
gzip ../SRR949205_2_val_2.fq

bismark --bowtie2 --phred33-quals --fastq -L 25 -N 1 --multicore 6 /home/shg047/db/hg19/bismark/ -1 SRR949205_1_val_1.fq.gz -2 SRR949205_2_val_2.fq.gz -o ./

find . -exec touch {} \;


Dear Sir/Madam,
This user use my photo and my name as its profile and post kinds of things full of defamation and insult and post my private things through google plus to public. It break the California laws and the profiles are with full of defamation and insult. Please remove all the photos and my private things. Please warn the user it is illegal. 

https://plus.google.com/100087098265453308274
 
https://plus.google.com/u/0/108664216580872773217
 
https://plus.google.com/collection/UpxPlB?hl=en-US
 
https://plus.google.com/collection/04KQlB?hl=en-US
 
https://plus.google.com/117041654657062516807?hl=en-US
 
https://plus.google.com/100087098265453308274?hl=en-US

If you have any questions, please call me.  my phone number is: 281-685-5882

By the way, This is not the first time she do that,  I have noticed google plus week. google plus delete her profiles. However, she create a fake profile to impersonate me again. Please remove the profile as soon as possible. 

Thanks. 




find . -maxdepth 5 -type f -exec touch {} +  &


cp /home/shg047/work/Minghua2016/fastq/N53_R1.fastq.gz /home/shg047/work/Minghua2016/fastq/N53_R2.fastq.gz /home/shg047/software/Bismark/test

mv /home/shg047/work/DennisLo2015/sortBam/HPLC.bam.bai ./HPLC.bismark.bam.bai
mv /home/shg047/work/DennisLo2015/sortBam/HPLC.bam ./HPLC.bismark.bam


### Produce Run Time
my $end_run = time();
my $run_time = $end_run - $start_run;
my $days  = int($run_time/(24*60*60));
my $hours = ($run_time/(60*60))%24;
my $mins  = ($run_time/60)%60;
my $secs  = $run_time%60;

warn "Bismark completed in ${days}d ${hours}h ${mins}m ${secs}s\n";
print REPORT "Bismark completed in ${days}d ${hours}h ${mins}m ${secs}s\n";

chr14:95233733-95237069

samtools view -b CTR107_trimmed.fq.gz_bismark_bt2.sort.bam chr14:95233733-95237069

samtools view -H file.bam > header.sam  # extract header only
samtools reheader header.sam file.unique.bam


/home/shg047/work/Alice/WNT5A/mouse/bam

samtools view -H SCNT.wnt5a_mouse.sort.bam > header.txt
samtools view SCNT.wnt5a_mouse.sort.bam | tail -n 300 > test.sam
samtools view -b -T ~/db/mm9/mm9.fa test.sam -o test.bam
samtools sort test.bam -o test.sort.bam
samtools index test.sort.bam
samtools view test.sort.bam

/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/test.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt 

/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/test.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt > ./SCNT.hapInfo.txt
/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/SCNT.wnt5a_mouse.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes MM9.CpG.positions.plus1.txt > SCNT.hapInfo.plus1.txt
/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/SCNT.wnt5a_mouse.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes MM9.CpG.positions.minus1.txt > SCNT.hapInfo.minus.txt

chr14:29321905

samtools tview SCNT.wnt5a_mouse.sort.bam ~/db/mm9/mm9.fa

Edit: An example script would be something of the form: 

samtools view -h foo.bam | awk '{if(found==0) {if($2=="SN:MT") {found=1; getline;}} print $0}' | samtools view -bSo foo.reformated.bam -


/home/shg047/oasis/db/mm9/MM9.CpG.positions.txt

samtools tview test.sort.bam ~/db/mm9/mm9.fa
samtools tview SCNT.wnt5a_mouse.sort.bam ~/db/mm9/mm9.fa


awk '{print $1,$2+1}' OFS="\t" /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt > MM9.CpG.positions.plus1.txt
awk '{print $1,$2-1}' OFS="\t" /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt > MM9.CpG.positions.minus1.txt



125194864
 
perl ~/bin/hapinfo2LDR2.pl ipsc.wnt5a chr14:29318659-29338701 < iPSC.hapInfo.txt
perl ~/bin/hapinfo2LDR2.pl scnt.wnt5a chr14:29318659-29338701 < SCNT.hapInfo.txt

bedtools coverage -d -abam SCNT.wnt5a_mouse.sort.bam -b wnt5a.bed > scnt.wnt5a.cov
bedtools coverage -d -abam iPSC.wnt5a_mouse.sort.bam -b wnt5a.bed > ipsc.wnt5a.cov

bismark_methylation_extractor --comprehensive -s --bedGraph SCNT.wnt5a_mouse.sort.bam -o /home/shg047/work/Alice/WNT5A/mouse/bam/rlt

PileOMeth extract -p 1 -q 1 --minDepth 1 --mergeContext ~/db/mm9/mm9.fa SCNT.wnt5a_mouse.sort.bam
 
 
 bam2fastq SCNT.wnt5a_mouse.sort.bam  
  bam2fastq iPSC.wnt5a_mouse.sort.bam 
 mv s_6_M_sequence.txt iPSC.fastq
 mv s_5_M_sequence.txt SCNT.fastq

 for i in iPSC.fastq SCNT.fastq
 do
 bismark  /home/shg047/oasis/db/mm9/ $i
 done
 
 
wc -l SCNT.wnt5a_mouse.sort_CpG.bedGraph

# sort the new bam files
for i in `ls *.WNT5A_human.bam`
do
# samtools sort $i -o $i.sort.bam
samtools index $i.sort.bam
done


# get the haplotype
perl ~/bin/bam2hapinfo

# calculate LDR


# get 



# get methylation level



cd /oasis/tscc/scratch/shg047/Alice/sortBam
/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/sortBam/human.bed /home/shg047/work/Alice/sortBam/SRR1286725_trimmed_bismark_bt2.sortc.bam.WNT5A_human.bam bismark /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt > ../hapinfo/SRR1286725.hapInfo.txt

/home/shg047/bin/mergedBam2hapInfoV2.pl /home/shg047/work/Alice/WNT5A/mouse/bam/wnt5a.bed /home/shg047/work/Alice/WNT5A/mouse/bam/SCNT.wnt5a_mouse.sort.bam bisReadMapper /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt > ../hapinfo/SCNT.hapInfo.txt

sh SRR1286295.job &
sh SRR1286296.job &
sh SRR1286297.job &
sh SRR1286298.job &
sh SRR1286299.job &

Last time you mentioned: We should focus on two aspects: (i) the set of marker regions that you identified for mapping tissue-of-origin and tumors; (ii) the method/algorithm that you take to use these markers for making prediction.


# assign to the nearest promoter
bedtools sort -i /media/Home_Raid1/shg047/NAS3/Alice/pathway/R1/ips_nt.R2.bed.txt > /media/Home_Raid1/shg047/NAS3/Alice/pathway/R1/ips_nt.R2.sort.bed.txt
bedtools intersect -wo -a /media/Home_Raid1/shg047/NAS3/Alice/pathway/R1/ips_nt.R2.sort.bed.txt -b /media/Home_Raid1/shg047/db/hg19/hg19_refGene.bed > /media/Home_Raid1/shg047/NAS3/Alice/pathway/R1/ips_nt.R2.sort.bed.reference.txt


wget -r --user=zhang --password='HDIIWpP' ftp://igm-storage1.ucsd.edu/161121_D00611_0401_BH3J7CBCXY_Zhang10X/

bismark --bowtie2 --phred33-quals --fastq -L 25 -N 1 --multicore 6 /home/shg047/db/hg19/bismark/ -1 ../fastq_trim/T9_R1_trimmed.fq.gz -2 ../fastq_trim/T9_R2_trimmed.fq.gz -o ../bam --basename T9_R1

gem3-mapper -p --bisulfite-mode -I GCA_000001405.15_GRCh38_no_alt_analysis_set_BS.gem -s 1 -p -M 4

for i in `ls *.fastq`
do
gzip $i
done
 
perl ~/bin/pairAutoBismark.pl ../PRJNA201480.txt 33 non


RASSF1	/home/shg047/oasis/monod/Fig1B/bam/RASSF1.61-WGBS.bam	/home/shg047/oasis/monod/Fig1B/bam/mhb.bed
SHOX2	/home/shg047/oasis/monod/Fig1B/bam/SHOX2.61-WGBS.bam	/home/shg047/oasis/monod/Fig1B/bam/mhb.bed

perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt submit bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt

samtools sort RASSF1.61-WGBS.bam -o RASSF1.61-WGBS.sort.bam
samtools sort SHOX2.61-WGBS.bam -o SHOX2.61-WGBS.sort.bam

samtools index RASSF1.61-WGBS.sort.bam
samtools index SHOX2.61-WGBS.sort.bam

mv RASSF1.61-WGBS.sort.bam RASSF1.61-WGBS.bam
mv SHOX2.61-WGBS.sort.bam SHOX2.61-WGBS.bam
mv RASSF1.61-WGBS.sort.bam.bai RASSF1.61-WGBS.bam.bai
mv SHOX2.61-WGBS.sort.bam.bai SHOX2.61-WGBS.bam.bai

chr3	50373531	50384063	chr3:50373531-50384063	
chr3	157811261	157826490	chr3:157811261-157826490

RASSF1: chr3:50373531-50384063 
SHOX2: chr3:157811261-157826490

#!/bin/csh
#PBS -q pdafm
#PBS -l nodes=1:ppn=1
#PBS -l walltime=8:00:00
#PBS -o GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.log
#PBS -e GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.err
#PBS -V
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd $PBS_O_WORKDIR
perl ~/bin/hapinfo2LDR2.pl SHOX2 chr3:157811261-157826490 < SHOX2.hapInfo.txt &
perl ~/bin/hapinfo2LDR2.pl RASSF1 chr3:50373531-50384063 < RASSF1.hapInfo.txt &

 
perl ~/bin/bam2hapInfo2PBS.pl mhb.bed submit bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt

/media/Home_Raid1/shg047/monod/Fig1b/bam/mhb.bed

scp shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/monod/Fig1b/bam/* ./


shg047@genomeMiner:~/work/db/mm9$ /media/Home_Raid1/shg047/work/db/mm9/mm9.refGene.bed^C



samtools view -bS -t ~/db/hg19/hg19.fa.fai RASSF1.61-WGBS.sam -o ~/RASSF1.61-WGBS.bam
samtools view -bS -t ~/db/hg19/hg19.fa.fai SHOX2.61-WGBS.sam -o ~/SHOX2.61-WGBS.bam

mv /media/Home_Raid1/shg047/RASSF1.61-WGBS.bam ./
mv /media/Home_Raid1/shg047/SHOX2.61-WGBS.bam ./



zcat SRR1648425_1.fastq.gz | head -n 40000 > ./test/SRR1648425_1.fastq.gz
zcat SRR1648425_2.fastq.gz | head -n 40000 > ./test/SRR1648425_2.fastq.gz

for i in {1..5}
do
perl HMHDection.pl 6-P-$i.hapInfo.txt 6-T-$i.hapInfo.txt excl2.txt 6-$i.txt
perl HMHDection.pl 7-P-$i.hapInfo.txt 7-T-$i.hapInfo.txt excl2.txt 7-$i.txt
done

cd /home/shg047/work/monod/hapinfo/december/mix
cat ../6-*.txt.txt > CRC-MHM.txt
cat ../7-*.txt.txt > LC-MHM.txt

perl -lane 'print if ! /LOC/' CRC-MHM.txt > CRC-MHM2.txt
perl -lane 'print if ! /LOC/' LC-MHM.txt > LC-MHM2.txt

awk 'NR>1 {print $1}' CRC-MHM2.txt > CRC-MHM3.txt
awk 'NR>1 {print $1}' LC-MHM2.txt > LC-MHM3.txt

perl -p -i -e "s/[:-]/\t/g" CRC-MHM3.txt
perl -p -i -e "s/[:-]/\t/g" LC-MHM3.txt
 
sort -u CRC-MHM3.txt > CRC-MHM4.txt
sort -u LC-MHM3.txt > LC-MHM4.txt

#!/bin/csh
#PBS -q pdafm
#PBS -l nodes=1:ppn=1
#PBS -l walltime=8:00:00
#PBS -o GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.log
#PBS -e GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.err
#PBS -V
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd $PBS_O_WORKDIR
# merge hapinfo for samples from same group
perl ~/bin/hapinfoMerge.pl NT.hapInfo.txt > NT.HapInfo.txt 
perl ~/bin/hapinfoMerge.pl CT.hapInfo.txt > CT.HapInfo.txt 
perl ~/bin/hapinfoMerge.pl H1.hapInfo.txt > H1.HapInfo.txt 
# then to do mhb calling
perl ~/bin/hapinfo2mhb.pl NT.hapInfo.txt.SumUniq 0.3 > NT.mhb &
perl ~/bin/hapinfo2mhb.pl H1.hapInfo.txt.SumUniq 0.3 > H1.mhb &
perl ~/bin/hapinfo2mhb.pl CT.hapInfo.txt.SumUniq 0.3 > CT.mhb &

perl ~/bin/hapinfo2mhb.pl W.hapInfo.txt.SumUniq 0.3 > W.mhb &

/media/Home_Raid1/shg047/work/monod/hapinfo/mix

cat H1.mhb NT.mhb > N.MHB
cat CT.mhb > C.MHB

bedtools sort -i N.MHB | wc -l
bedtools sort -i C.MHB | wc -l
bedtools intersect -wa -v -a C.MHB -b N.MHB | sort -u | wc -l 

bedtools intersect -wa -a C.MHB -b N.MHB | sort -u | wc -l 
bedtools intersect -wa -v -a C.MHB -b N.MHB | sort -u | wc -l 
bedtools intersect -wa -v -a N.MHB -b C.MHB | sort -u | wc -l 

cd /home/shg047/work/monod/hapinfo/June
load("monod.mhl.22July.RData")

NT.hapInfo.txt
-rwxrwxrwx 1 dinh sambashare  629 Nov 23 13:49 hapinfo2mhb.job
-rwxrwxrwx 1 dinh sambashare 255M Nov 23 13:49 CT.hapInfo.txt
-rwxrwxrwx 1 dinh sambashare    0 Nov 23 13:49 NT.mhb
-rwxrwxrwx 1 dinh sambashare 508M Nov 23 13:49 H1.hapInfo.txt

cd /home/shg047/work/Alice/sortBam
perl ~/bin/SaminfoPre4hapinfo.pl > sampleconfig.txt
perl ~/bin/bam2hapInfo2PBS.pl sampleconfig.txt non bismark /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt


ls SRR1286395*
/home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed


#!/usr/bin/sh
START=$(date +%s)
END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
echo "It takes DIFF=$DIFF seconds to complete this task..."

ls *gz | awk -F. '{print $1}'
 
 diff <(grep "run complete" *.err | awk -F: '{print $1}'|sort) <(grep kill *.err | awk -F: '{print $1}' | sort)
 
grep "run complete" *.err | awk -F: '{print $1}'|sort > a
ls ../bam/*bam | grep -v temp | awk -F[/_] '{print $3".err"}' | sort | wc -l> b 
diff <(grep "run complete" *.err | awk -F: '{print $1}'|sort) <(ls ../bam/*bam | grep -v temp | awk -F[/_] '{print $3".err"}' )
paste a b > c


SRR1286589.pbs



xx_trimmed_bismark_bt2.bam

# bash date
START=$(date +%s)
END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
echo "It takes DIFF=$DIFF seconds to complete this task..."﻿

# perl time
my $start_run = time();
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";

Note: setting of bismark alignment in PBS system (DNA methylation, BS-seq, RRBS)

multicore=1, ppn=8, memory=3.2G*6
multicore=2, ppn=16, memory=3.2G*12

note: '-p 1' will already use 
* 4 threads/cores for Bowtie2 (3.2G/per) plus 1 additional core for Bismark itself
* 2 threads/cores for Perl (3.2G/per)
* 2 threads/cores for SAMTOOLS

note: '-p 2' will already use 
* 8 threads/cores for Bowtie2 (3.2G/per) plus 1 additional core for Bismark itself
* 4 threads/cores for Perl (3.2G/per)
* 4 threads/cores for SAMTOOLS

62067 shg047    20   0 3353m 3.1g  812 R 100.0  5.0  32:53.22 perl
61902 shg047    20   0 3353m 3.1g 2168 R 100.0  5.0  33:48.64 perl
62048 shg047    20   0 3353m 3.1g 2168 R 100.0  5.0  33:29.31 perl
62126 shg047    20   0 3353m 3.1g  812 R 99.7  5.0  32:47.26 perl
62155 shg047    20   0 3357m 3.2g 1968 S 59.3  5.0  18:00.17 bowtie2-align-s
62154 shg047    20   0 3357m 3.2g 1968 S 58.9  5.0  17:55.70 bowtie2-align-s
62169 shg047    20   0 3357m 3.2g 1968 S 58.6  5.0  17:46.12 bowtie2-align-s
62146 shg047    20   0 3357m 3.2g 1968 S 58.3  5.0  17:55.60 bowtie2-align-s
62167 shg047    20   0 3357m 3.2g 1968 S 57.6  5.0  17:45.36 bowtie2-align-s
62171 shg047    20   0 3357m 3.2g 1968 S 57.6  5.0  17:50.46 bowtie2-align-s
62149 shg047    20   0 3357m 3.2g 1968 S 55.6  5.0  18:00.27 bowtie2-align-s
62168 shg047    20   0 3357m 3.2g 1968 S 55.3  5.0  17:49.37 bowtie2-align-s
62178 shg047    20   0 22884 5388  780 S 13.9  0.0   4:07.75 samtools
62183 shg047    20   0 22884 5388  780 S 13.6  0.0   4:08.99 samtools
62180 shg047    20   0 22884 5388  780 S 13.2  0.0   4:08.05 samtools
62184 shg047    20   0 22884 5388  780 S 13.2  0.0   4:08.89 samtools



/home/shg047/oasis/db/hg19/align/bismark

/media/Home_Raid1/shg047/git/Bismark/bismark ~/db/hg19/align/bismark/ --multicore=4 --non_directional -1 N36_R1.fastq -2 N36_R2.fastq

perl /media/Home_Raid1/shg047/git/Bismark/


# How to install git with adminstration 
sudo apt-get update
sudo apt-get install git
git config --global user.name "Shicheng-Guo"
git config --global user.email "shicheng.guo@hotmail.com"
git config --list
git commit --amend --reset-author

# How to install git without admin (TSCC in UCSD)
cd software
wget https://www.kernel.org/pub/software/scm/git/git-2.10.2.tar.xz
tar xf xzvf git-2.10.2.tar.xz
cd git-2.10.2
./configure --prefix="$HOME/git"
make install
git config --list

# Go to home and bulid the databse for git
cd
mkidr git



	
cd 
scp /media/Home_Raid1/shg047/NAS3/Alice/human/fastq
scp ../human/SRR1286624.fastq.gz.bismark.pbs shg047@genome-miner.ucsd.edu:/media/Home_Raid1/shg047/NAS3/Alice/human/fastq

ls *fastq | echo

ls *fastq | paste - -




wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/pysam/pysam-0.6.tar.gz

cd /media/Home_Raid1/shg047/NAS3/Minghua2016/fastq
bedtools unionbedg -i *bedgraph -filler NA -header -names CT1.bam.bedgraph CT2.bam.bedgraph CT3.bam.bedgraph CT.bam.bedgraph CU1.bam.bedgraph CU2.bam.bedgraph CU3.bam.bedgraph CU.bam.bedgraph ET1.bam.bedgraph ET2.bam.bedgraph ET3.bam.bedgraph ET.bam.bedgraph EU1.bam.bedgraph EU2.bam.bedgraph EU3.bam.bedgraph EU.bam.bedgraph HE.bam.bedgraph N100.bam.bedgraph N101.bam.bedgraph N102.bam.bedgraph N103.bam.bedgraph N10.bam.bedgraph N11.bam.bedgraph N12.bam.bedgraph N13.bam.bedgraph N14.bam.bedgraph N15.bam.bedgraph N16.bam.bedgraph N17.bam.bedgraph N18.bam.bedgraph N19.bam.bedgraph N1.bam.bedgraph N20.bam.bedgraph N21.bam.bedgraph N22.bam.bedgraph N23.bam.bedgraph N24.bam.bedgraph N25.bam.bedgraph N26.bam.bedgraph N27.bam.bedgraph N28.bam.bedgraph N29.bam.bedgraph N2.bam.bedgraph N30.bam.bedgraph N31.bam.bedgraph N32.bam.bedgraph N33.bam.bedgraph N34.bam.bedgraph N35.bam.bedgraph N36.bam.bedgraph N37.bam.bedgraph N38.bam.bedgraph N39.bam.bedgraph N3.bam.bedgraph N40.bam.bedgraph N41.bam.bedgraph N42.bam.bedgraph N43.bam.bedgraph N44.bam.bedgraph N45.bam.bedgraph N46.bam.bedgraph N47.bam.bedgraph N48.bam.bedgraph N49.bam.bedgraph N4.bam.bedgraph N50.bam.bedgraph N51.bam.bedgraph N52.bam.bedgraph N53.bam.bedgraph N54.bam.bedgraph N55.bam.bedgraph N56.bam.bedgraph N57.bam.bedgraph N58.bam.bedgraph N59.bam.bedgraph N5.bam.bedgraph N60.bam.bedgraph N61.bam.bedgraph N62.bam.bedgraph N63.bam.bedgraph N64.bam.bedgraph N65.bam.bedgraph N66.bam.bedgraph N67.bam.bedgraph N68.bam.bedgraph N69.bam.bedgraph N6.bam.bedgraph N70.bam.bedgraph N71.bam.bedgraph N72.bam.bedgraph N73.bam.bedgraph N74.bam.bedgraph N75.bam.bedgraph N76.bam.bedgraph N77.bam.bedgraph N78.bam.bedgraph N79.bam.bedgraph N7.bam.bedgraph N80.bam.bedgraph N81.bam.bedgraph N82.bam.bedgraph N83.bam.bedgraph N84.bam.bedgraph N85.bam.bedgraph N86.bam.bedgraph N87.bam.bedgraph N88.bam.bedgraph N89.bam.bedgraph N8.bam.bedgraph N90.bam.bedgraph N91.bam.bedgraph N92.bam.bedgraph N93.bam.bedgraph N94.bam.bedgraph N95.bam.bedgraph N96.bam.bedgraph N97.bam.bedgraph N98.bam.bedgraph N99.bam.bedgraph N9.bam.bedgraph T100.bam.bedgraph T101.bam.bedgraph T102.bam.bedgraph T103.bam.bedgraph T10.bam.bedgraph T11.bam.bedgraph T12.bam.bedgraph T13.bam.bedgraph T14.bam.bedgraph T15.bam.bedgraph T16.bam.bedgraph T17.bam.bedgraph T18.bam.bedgraph T19.bam.bedgraph T1.bam.bedgraph T20.bam.bedgraph T21.bam.bedgraph T22.bam.bedgraph T23.bam.bedgraph T24.bam.bedgraph T25.bam.bedgraph T26.bam.bedgraph T27.bam.bedgraph T28.bam.bedgraph T29.bam.bedgraph T2.bam.bedgraph T30.bam.bedgraph T31.bam.bedgraph T32.bam.bedgraph T33.bam.bedgraph T34.bam.bedgraph T35.bam.bedgraph T36.bam.bedgraph T37.bam.bedgraph T38.bam.bedgraph T39.bam.bedgraph T3.bam.bedgraph T40.bam.bedgraph T41.bam.bedgraph T42.bam.bedgraph T43.bam.bedgraph T44.bam.bedgraph T45.bam.bedgraph T46.bam.bedgraph T47.bam.bedgraph T48.bam.bedgraph T49.bam.bedgraph T4.bam.bedgraph T50.bam.bedgraph T51.bam.bedgraph T52.bam.bedgraph T53.bam.bedgraph T54.bam.bedgraph T55.bam.bedgraph T56.bam.bedgraph T57.bam.bedgraph T58.bam.bedgraph T59.bam.bedgraph T5.bam.bedgraph T60.bam.bedgraph T61.bam.bedgraph T62.bam.bedgraph T63.bam.bedgraph T64.bam.bedgraph T65.bam.bedgraph T66.bam.bedgraph T67.bam.bedgraph T68.bam.bedgraph T69.bam.bedgraph T70.bam.bedgraph T71.bam.bedgraph T72.bam.bedgraph T73.bam.bedgraph T74.bam.bedgraph T75.bam.bedgraph T76.bam.bedgraph T77.bam.bedgraph T78.bam.bedgraph T79.bam.bedgraph T7.bam.bedgraph T80.bam.bedgraph T81.bam.bedgraph T82.bam.bedgraph T83.bam.bedgraph T84.bam.bedgraph T85.bam.bedgraph T86.bam.bedgraph T87.bam.bedgraph T88.bam.bedgraph T89.bam.bedgraph T8.bam.bedgraph T90.bam.bedgraph T91.bam.bedgraph T92.bam.bedgraph T93.bam.bedgraph T94.bam.bedgraph T95.bam.bedgraph T96.bam.bedgraph T97.bam.bedgraph T98.bam.bedgraph T99.bam.bedgraph T9.bam.bedgraph > esca.bg
python ~/software/BSseeker2/bs_seeker2-align.py -1 N36_R1.fastq -2 N36_R2.fastq --aligner=bowtie2 -o ./96.bam -f bam -g /media/Home_Raid1/shg047/NAS3/Minghua2016/fa/target.fa
python ~/software/BSseeker2/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i 96.bam -o 96 --txt --db /media/Home_Raid1/shg047/software/BSseeker2/bs_utils/reference_genomes/target.fa_bowtie2


bedtools interact /home/shg047/db/hg19/lncRNA.hg19.bed
132.239.25.238

awk '{print $1,$3,$4,$10}' -OFS="\t" gencode.v19.long_noncoding_RNAs.gtf | head

bismark -q --phred33-quals --multicore=1 -n 1 -l 20 --non_directional /media/Home_Raid1/shg047/db/aligndb/hg19/bismark N36.R1.fastq -o ../bam/

chr1	35037044	35037796	chr1:35037044-35037796

Sylvia Schmalz 

for i in `ls *bam`
do

python /media/Home_Raid1/shg047/software/BSseeker2/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i $i -o $i --txt --db /home/shg047/software/BSseeker2/bs_utils/reference_genomes/hg19.fa_bowtie2
done




ls -larth SRR1286295.fastq.gz
ls -larth SRR1286296.fastq.gz 
ls -larth SRR1286297.fastq.gz
ls -larth SRR1286298.fastq.gz
ls -larth SRR1286299.fastq.gz


ENSG00000270170.1

GSM978967
cat cerebellum_MethylC-seq_chr*.BED > GSM978968.Brain.BED

for i in {1..11}
do
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer$i/tracks_hg19/Human_PancreaticCancer$i.meth.bw
done
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas1/tracks_hg19/Human_NormalPancreas1.meth.bw
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas2/tracks_hg19/Human_NormalPancreas2.meth.bw
wget http://smithlab.usc.edu/methbase/data/Gao-Human-2015/Human_BloodHealthy/tracks_hg19/Human_BloodHealthy.meth.bw

Step 1: install GSL
wget http://mirror.keystealth.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd /home/shg047/software/gsl-2.1
./configure --prefix=/media/Home_Raid1/shg047/software/gsl-2.2.1
make
make install
export CPATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/include
export LIBRARY_PATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/lib

Step2： install methpipe
wget http://smithlabresearch.org/downloads/methpipe-3.4.2.tar.bz2
tar xjvf methpipe-3.4.2.tar.bz2
make
﻿export PATH=$PATH:/media/Home_Raid1/shg047/software/methpipe-3.4.2/bin

Step 3: install rmap
wget http://smithlabresearch.org/downloads/rmap-2.1.tar.bz2
tar xjvf rmap-2.1.tar.bz2
cd rmap-2.1/
make
sudo make install
PATH=$PATH:/media/Home_Raid1/shg047/software/rmap-2.1/bin

Step 4. Install walt
git clone https://github.com/smithlabcode/walt.git
cd walt
make

Step 5: Download Genome Reference - hg19
cd /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

Step 6: 
makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex






/media/Home_Raid1/shg047/NAS3/tcga/chol/lncrna/bam
cd /home/shg047/oasis/Minghua2016/fastq
bismark -q --phred33-quals -n 1 -l 20 --non_directional /home/shg047/oasis/db/hg19 -1 CT1_R1.fastq -2 CT1_R2.fastq -o ../

bismark -q --phred33-quals --multicore=8 -n 1 -l 20 --non_directional /media/Home_Raid1/shg047/db/aligndb/hg19/bismark -1 T96.R1.fastq -2 T96.R2.fastq -o ./
bismark -q --phred33-quals --multicore=8 -n 1 -l 20 --non_directional /media/Home_Raid1/shg047/db/aligndb/hg19/bismark -1 N36.R1.fastq -2 N36.R2.fastq -o ./



wget http://pellegrini.mcdb.ucla.edu/BS_Seeker2/sup_data/DS2_PE_simu_perfect_end1.fa.gz
wget http://pellegrini.mcdb.ucla.edu/BS_Seeker2/sup_data/DS2_PE_simu_perfect_end2.fa.gz

pip install 'pysam=6,<7'

cd /home/shg047/oasis/Minghua2016/fastq
python /home/shg047/software/BSseeker2-master/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i N46.bam -o N46.mr --txt --db /home/shg047/software/BSseeker2/bs_utils/reference_genomes/hg19.fa_bowtie2
bs_seeker2


cd /home/shg047/db/hg19/
REF="/home/shg047/db/hg19/hg19.fa";
python /home/shg047/software/bwa-meth/bwameth.py index $REF
REF="/home/shg047/db/hg19/hg19.fa";
python /home/shg047/software/bwa-meth/bwameth.py --reference $REF N28_R1.fastq N28_R2.fastq | samtools view -b - > N28.bwameth.bam
samtools view -cf2 N28.bwameth.bam
samtools view -cF4 N28.bwameth.bam
samtools flagstat N28.bwameth.bam

/home/shg047/software/Python-2.7.4/python /home/shg047/software/BSseeker2_v2.0.9/bs_seeker2-call_methylation.py -x -r 5 --rm-CCGG -i N46.bam -o N46.mr --txt --db /home/shg047/software/BSseeker2/bs_utils/reference_genomes/hg19.fa_bowtie2
bs_seeker2


	 samtools index T19.bam
	 samtools sort T19.bam -o T19.sort
	 

wget ftp://ftp.gnu.org/pub/gnu/emacs/emacs-25.1.tar.xzvf
tar xf emacs-25.1.tar.xzvf
./configure
make
export PATH=$PAHT:

cd /home/shg047/oasis/Minghua2016/fastq

#!/usr/bin/sh
mkdir $HOME/db
cd $HOME/db
for i in hg18 hg19 hg38 mm9 mm10 
do
[! -d "$i"] && mkdir $i 
cd $i
mkdir "fa"
cd "fa"
wget -r -np http://hgdownload.soe.ucsc.edu/goldenPath/$i/bigZips/$i.chromFa.tar.gz  
wget -r -np http://hgdownload.soe.ucsc.edu/goldenPath/$i/bigZips/chromFa.tar.gz     
wget -r -np http://hgdownload.soe.ucsc.edu/goldenPath/$i/bigZips/$i.chrom.sizes
tar xzvf *.gz
rm *.gz
cat *fa ../$i.fa
cd ../../
done

wget -r http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/$i.chromFa.tar.gz  
wget -r http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz  
tar xzvf *.gz
rm *.gz
cat *fa ../$i.fa


python ~/software/BSseeker2/bs_seeker2-build.py -f /home/shg047/db/hg19/hg19.fa --aligner=bowtie2
/home/shg047/software/Python-2.7.12/python /home/shg047/software/BSseeker2/bs_seeker2-align.py -1 CT1_R1.fastq -2 CT1_R2.fastq -g /home/shg047/oasis/db/hg19/hg19.fa -t Y
cd 

python ~/software/BSseeker2/bs_seeker2-build.py -f /home/puweilin/Software/bowtie2-2.2.9/ChrM/Homo_sapiens.GRCh38.dna.chromosome.MT.fa --aligner=bowtie2

python ~/software/BSseeker2/bs_seeker2-align.py -1 T86_R1.fastq -2 T86_R2.fastq --aligner=bowtie2 -o ./T86.bam -f bam -g /home/shg047/db/hg19/hg19.fa

# genome-miner
/home/shg047/software/Python-2.7.4/python ~/software/BSseeker2/bs_seeker2-build.py -f /home/shg047/db/hg19/hg19.fa --aligner=bowtie2


Step 1: install GSL
wget http://mirror.keystealth.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd /home/shg047/software/gsl-2.1
./configure --prefix=/media/Home_Raid1/shg047/software/gsl-2.2.1
make
make install
export CPATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/include
export LIBRARY_PATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/lib

Step2： install methpipe
wget http://smithlabresearch.org/downloads/methpipe-3.4.2.tar.bz2
tar xjvf methpipe-3.4.2.tar.bz2
make
﻿export PATH=$PATH:/media/Home_Raid1/shg047/software/methpipe-3.4.2/bin

Step 3: install rmap
wget http://smithlabresearch.org/downloads/rmap-2.1.tar.bz2
tar xjvf rmap-2.1.tar.bz2
cd rmap-2.1/
make
sudo make install
PATH=$PATH:/media/Home_Raid1/shg047/software/rmap-2.1/bin

Step 4. Install walt
git clone https://github.com/smithlabcode/walt.git
cd walt
make

Step 5: Download Genome Reference - hg19
cd /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

Step 6: 
makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex





bsrate -c ~/oasis/db/hg19/hg19.fa -o T84.mr.bsrate T84.mr

# Line-1 Methylation Primer Set-1
5′-GGACGTATTTGGAAAATCGGG-3′ 
5′-AATCTCGCGATACGCCGTT-3′ 
5′-TCGAATATTGCGTTTTCGGATCGGTTT-3′ 

# Line-1 Methylation Primer Set-2
5′-TTGAGTTGTGGTGGGTTTTATTTAG-3′ 
5′-TCATCTCACTAAAAAATACCAAACA-3′.

# Line-1 Methylation Primer Set-3 
5’-TTGGTTAGGTGTGGGATATAGTT-3’
5’-CAAAAAATCAAAAAATTCCCTTTCC-3’

mkdir pip
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py --prefix="./"
pip install --user --upgrade cutadapt

~/.local/bin/cutadapt --help
cp ~/.local/bin/cutadapt ~/bin/
cp ~/.local/bin/cutadapt /media/Home_Raid1/shg047/software/trim_galore_zip

git clone https://github.com/BSSeeker/BSseeker2.git

git clone https://github.com/marcelm/cutadapt.git
bismark -q --phred33-quals -n 1 -l 40 --non_directional ~/NAS3/
bismark -q --phred33-quals --bowtie2 -N 1 -L 30 ~/NAS3/db/hg19/ -1 T84_R1.fastq -2 T84_R2.fastq
bismark -q --phred33-quals -n 1 -l 30 ~/NAS3/db/hg19/ -1 T84_R1_val_1.fq -2 T84_R2_val_2.fq
bismark --bowtie2 --phred33-quals --RRBS -- --fastq -L 30 -N 1 --multicore 4 ~/NAS3/db/hg19/ -1 T84_R1_val_1.fq -2 T84_R2_val_2.fq  -o ./

trim_galore --phred33 --paired --fastqc --non_directional --rrbs --illumina -o ./ T84_R1.fastq T84_R2.fastq
trim_galore --paired -a GATCGGAAGAGCA -a2 GCTCTTCCGATCT --retain_unpaired  --trim1  T21.5.read1.fq.head T21.5.read2.fq.head  

chr3:142839987
chr17:59534637
chr2:74742664
chr19:58446159-58447359
19	58446371
chr12:95942887

walt -i ~/NAS3/db/hg19/chrosome/methpipe.hg19.dbindex -1 T84_R1.fastq -2 T84_R2.fastq -o T84.mr


for i in `ls *gz`
do
gunzip -f $i
done

bismark --bowtie2 --phred33-quals --fastq -L 30 -N 1 --multicore 4 ~/NAS3/db/hg19/ -1 T83_R1.fastq.gz -2 T83_R2.fastq.gz  -o ../bam
chr10:118810290
chr13:103498159
chr19:37341777
chr19:12163473
chr10:123922851
chr12:128751882

makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex

# 2016-11-03
cat cerebellum_MethylC-seq_chr*.BED > GSM978968.Brain.BED

for i in {1..11}
do
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer$i/tracks_hg19/Human_PancreaticCancer$i.meth.bw
done
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas1/tracks_hg19/Human_NormalPancreas1.meth.bw
wget http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_NormalPancreas2/tracks_hg19/Human_NormalPancreas2.meth.bw
wget http://smithlab.usc.edu/methbase/data/Gao-Human-2015/Human_BloodHealthy/tracks_hg19/Human_BloodHealthy.meth.bw

# 2016-11-03
Step 1: install GSL
wget http://mirror.keystealth.org/gnu/gsl/gsl-latest.tar.gz
tar xzvf gsl-latest.tar.gz
cd /home/shg047/software/gsl-2.1
./configure --prefix=/media/Home_Raid1/shg047/software/gsl-2.2.1
make
make install
export CPATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/include
export LIBRARY_PATH=/media/Home_Raid1/shg047/software/gsl-2.2.1/lib

Step2： install methpipe
wget http://smithlabresearch.org/downloads/methpipe-3.4.2.tar.bz2
tar xjvf methpipe-3.4.2.tar.bz2
make
﻿export PATH=$PATH:/media/Home_Raid1/shg047/software/methpipe-3.4.2/bin

Step 3: install rmap
wget http://smithlabresearch.org/downloads/rmap-2.1.tar.bz2
tar xjvf rmap-2.1.tar.bz2
cd rmap-2.1/
make
sudo make install
PATH=$PATH:/media/Home_Raid1/shg047/software/rmap-2.1/bin

Step 4. Install walt
git clone https://github.com/smithlabcode/walt.git
cd walt
make

Step 5: Download Genome Reference - hg19
cd /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

Step 6: 
makedb -c /media/Home_Raid1/shg047/NAS3/db/hg19/chrosome -o methpipe.hg19.dbindex



for i in `ls *bw`
do
bigWigAverageOverBed $i ES.bed $i.es.out
done

# download histone modification from roadmap/encode
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29611/
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56712
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56712/suppl/GSE56712_RAW.tar
tar xvf GSE56712_RAW.tar
rm GSE56712_RAW.tar
gunzip *.gz
rm *bam
rm *bigWig

# In terms of Super-Enhancer(H3k27ac)
cat *H3k27ac*broadPeak > Enhancer.txt
awk '{if ($7>5 && $5>100) print $1,$2,$3}' OFS="\t" Enhancer.txt > Enhancer.bed
sort -k1,1 -k2,2n Enhancer.bed > EnhancerSort.bed
bedtools merge -i EnhancerSort.bed > Human-Enhancer-hg19-H3k27ac-PMID24119843.bed
mv Human-Enhancer-hg19-PMID24119843.bed Human-Enhancer-hg19-H3k27ac-PMID24119843.bed
wc -l Human-Enhancer-hg19-H3k27ac-PMID24119843.bed

perl bed3Tobed4.pl Human-hg19-H3k27ac-PMID24119843.bed > Human-hg19-H3k27ac-PMID24119843.bed4

bedtools intersect -u -wa -a ~/NAS3/db/GPL13534.map -b superEnhancerCluster.bed | wc -l
bedtools intersect -u -wa -a superEnhancerCluster.bed -b ~/NAS3/db/GPL13534.map | wc -l
bedtools intersect -u -wa -a superEnhancerCluster.bed -b ~/NAS3/db/hg19/CpGI.hg19.bed | wc -l
bedtools intersect -u -wa -a superEnhancerCluster.bed -b ~/NAS3/db/hg19/CpG.Shore.hg19.bed | wc -l
bedtools intersect -u -wa -a superEnhancerCluster.bed -b ~/NAS3/db/hg19/CpG.Shelf.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shore.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shelf.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpGI.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shelf.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpGI.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpG.Shore.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpG.Shelf.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/Enhancers.Fantom.hg19.bed -b superEnhancer.bed | wc -l

# In terms of Promoter(H3k4me3)
zcat *H3k4me3*broadPeak.gz > promter.txt
awk '{if ($7>5) print $1,$2,$3}' OFS="\t" promter.txt > promter.bed
sort -k1,1 -k2,2n promter.bed > promterSort.bed
bedtools merge -i promterSort.bed > promterSortMerge.bed
wc -l promterSortMerge.bed
bedtools intersect -u -wa -a ~/NAS3/db/GPL13534.map -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpGI.hg19.bed -b promterSortMerge.bed| wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shore.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shelf.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpGI.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/CpG.Shelf.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b promterSortMerge.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpGI.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpG.Shore.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.all.hg19.bed -b ~/NAS3/db/hg19/CpG.Shelf.hg19.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/enhancer.12ct.encode.hg19.bed -b superEnhancer.bed | wc -l
bedtools intersect -u -wa -a ~/NAS3/db/hg19/Enhancers.Fantom.hg19.bed -b superEnhancer.bed | wc -l

TCGA-ESCA-PBMC.bed
bedtools intersect -u -wa -a TCGA-ESCA-PBMC.bed -b ~/NAS3/db/hg19/ 




ls -larth *H3k27ac*broadPeak | wc -l
ls -larth *H3k27me3*broadPeak | wc -l
ls -larth *H3k79me2*broadPeak | wc -l
ls -larth *H3k9me3*broadPeak | wc -l
ls -larth *H3k4me2*broadPeak | wc -l

cd 
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56712/suppl/GSE56712_RAW.tar
tar xvf GSE56712_RAW.tar
rm GSE56712_RAW.tar
rm *bam
rm *bigWig
gunzip *.gz
for i in H3k4me1 H3k4me2 H3k4me3 H3k9ac H3k9me1 H3k9me3 H3k27ac H3k27me3 H3k36me3 H3k79me2 H4k20me1 
do
cat *$i*broadPeak > promter.txt
awk '{if ($7>5) print $1,$2,$3}' OFS="\t" promter.txt > promter.bed
sort -k1,1 -k2,2n promter.bed > promterSort.bed
bedtools merge -i promterSort.bed > Human-hg19-$i-PMID24119843.bed
done
rm promter.txt
rm promter.bed
rm promterSort.bed


TCGA-ESCA-PBMC.bed

cd /media/Home_Raid1/shg047/NAS3/HM450/TCGA/esca



#!/bin/csh
#PBS -q hotel
#PBS -l nodes=1:ppn=1
#PBS -l walltime=168:00:00
#PBS -o GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.log
#PBS -e GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.err
#PBS -V
#PBS -M shicheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd /oasis/tscc/scratch/shg047/Roadmap/wig
bigWigCorrelate bwList.txt > bwListCorrelationResult.txt


for i in chr10:79857968-79858049 chr10:79858088-79858205 chr10:79858211-79858298 chr10:79858319-79858347 chr6:122658641-122658674 chr6:122659252-122659297
do
perl ~/bin/hapinfo2LDR2.pl $i.tissue.R2.txt $i < tissue.hapinfo.txt 
perl ~/bin/hapinfo2LDR2.pl $i.ips.R2.txt $i < ipsNT.hapinfo.txt
perl ~/bin/hapinfo2LDR2.pl $i.esc.R2.txt $i < ESC.hapinfo.txt
done

file<-list.files(pattern="*.rsq$")
for(i in file){
print(i)
M<-read.table(i,head=T,row.names=1,as.is=T)
library("grDevices")
col=colorRampPalette(c("white", "red"))(10) 
M[lower.tri(M)] <- NA
pdf(paste(i,"pdf",sep="."))
image(data.matrix(M),col = col,frame=F,xaxt="n",yaxt="n")
dev.off()
}

#!/bin/csh
#PBS -q hotel
#PBS -l nodes=1:ppn=1
#PBS -l walltime=9:00:00
#PBS -o GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.log
#PBS -e GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.err
#PBS -V
#PBS -M shicheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd /oasis/tscc/scratch/shg047/Roadmap/wig
wigToBigWig GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz ~/oasis/db/hg19/hg19.chrom.sizes GSM1010980_UCSD.Ovary.Bisulfite-Seq.STL002.wig.gz.bw


/home/shg047/oasis/Holger2016/bedGraph

cd ~/bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed  ./
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat ./
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCluster ./
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCorrelate ./
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigInfo ./
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigMerge ./
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary ./
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph ./
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig ./


for i in `ls GSM*`
do
perl GSE52270.pl $i > $i.bedgraph &
done

for i in `ls *bedgraph`
do
sort -k1,1 -k2,2n $i > $i.sort &
done

for i in `ls *bedgraph`
do
bedGraphToBigWig $i.sort ~/oasis/db/hg19/hg19.chrom.sizes $i.sort.hg19.bw &
done

for i in `ls *.wig.gz`
do
wigToBigWig $i ~/oasis/db/hg19/hg19.chrom.sizes $i.bw &
done





for i in `ls *hapinfo.txt`
do
perl ~/bin/hapinfo2BlocAvgR2.pl $i > R2.$i
done
perl ../

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl tissue.hapinfo.txt > tissue.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl tissue.hapinfo.txt.$i > R2.tissue.hapinfo.txt.$i
done

for i in {1..100}
do
perl ~/bin/randomSampleFromHaploInfo.pl tissue.hapinfo.txt > tissue.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl tissue.hapinfo.txt.$i > R2.tissue.hapinfo.txt.$i
perl ~/bin/randomSampleFromHaploInfo.pl ipsNT.hapinfo.txt > ipsNT.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ipsNT.hapinfo.txt.$i > R2.ipsNT.hapinfo.txt.$i
perl ~/bin/randomSampleFromHaploInfo.pl ESC.hapinfo.txt > ESC.hapinfo.txt.$i
perl ~/bin/hapinfo2BlocAvgR2.pl ESC.hapinfo.txt.$i > R2.ESC.hapinfo.txt.$i
done

cp /oasis/tscc/scratch/zhl002/WGBS_mouse/HapInfo/ESC.hapinfo.txt.SumUniq  /home/shg047/oasis/mouse/hapinfo/group/ESC.hapinfo.txt 
cp /oasis/tscc/scratch/zhl002/WGBS_mouse/HapInfo/ipsNT.hapinfo.txt.SumUniq  /home/shg047/oasis/mouse/hapinfo/group/ipsNT.hapinfo.txt 
cp /oasis/tscc/scratch/zhl002/WGBS_mouse/HapInfo/tissue.hapinfo.txt.SumUniq  /home/shg047/oasis/mouse/hapinfo/group/tissue.hapinfo.txt

perl 


cd /home/shg047/oasis/mouse/sortBam
perl ../saminfoPre4bam2hapinfo.pl > Samconfig2.txt
perl /home/shg047/bin/bam2hapInfo2PBS.pl Samconfig2.txt submit bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt

for i in `ls *fastq`
do
zip $i.gz $i
done


perl ../saminfoPre4bam2hapinfo.pl > SamConfig.txt

SRX209451
SRX271137
SRX271142


 #!/bin/csh
 #PBS -N mf2bedGraph
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
 perl ~/bin/hapinfo2LDR2.pl  rlt chr10:110930425-110930618 Mouse.Merge.Hapinfo.txt.SumUniq
 

 #!/bin/csh
 #PBS -N mf2bedGraph
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
 perl ~/bin/hapinfo2mhb.pl Mouse.Merge.Hapinfo.txt.SumUniq 0.5 > MHB.RD90up80.r0.5-2.bed

 
 cd /home/shg047/oasis/mouse/hapinfo
perl /home/shg047/oasis/mouse/mergeHapinfo/hapinfoMerge.pl Mouse.Merge.Hapinfo.txt Mouse.Merge.Hapinfo.Merge.txt

 #!/bin/csh
 #PBS -N mf2bedGraph
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
 perl /home/shg047/oasis/mouse/mergeHapinfo/hapinfoMerge.pl Mouse.Merge.Hapinfo.txt Mouse.Merge.Hapinfo.Merge.txt
 
 #!/bin/csh
 #PBS -N hapinfo2mhl
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
perl ~/bin/MethylFreq2wig.pl 
 
 
cd /oasis/tscc/scratch/zhl002/WGBS_mouse/bam
chr1:3000017-3000744
   
   view -b Indx01.merged.bam.sorted.bam chr1:3000017-3000744 > xx.bam
   
   
qsub SRX080192.job
qsub SRX1091397.job
qsub SRX209448.job
qsub SRX209453.job
qsub SRX209458.job
qsub SRX271136.job



/home/shg047/oasis/mouse/mergeHapinfo/hapinfo/MHB.Mouse.R2-0.5.bed



for i in `ls *job`
do
qsub $i
done

rm SRX1019864.job
rm SRX1019865.job
rm SRX1019866.job
rm SRX1019867.job

for i in {4..7}
do
qsub SRX101986$i.job
done




samtoos view SRR630257_1_trimmed_bismark_bt2.bam
SRR630258_1_trimmed_bismark_bt2.bam

SRR630257.12464425_HWI-ST1113_0085:5:1308:18692:107279_length=100       16      chr4    2999991 8       95M     *       0       0       TCACCAACAAAAATTCACTCGAACAACTAAATTCTTCTTTTTTTTTTTCCATT
SRR630257.5426488_HWI-ST1113_0085:5:1203:8069:184455_length=100 0       chr4    3001064 42      95M     *       0       0       GATTTTTAGAGTGGTTGTATAAGTTTAGAATTTTATTAATAATGGAAGAGTGTTTTTTTTT

SRR630258.15688879_HWI-ST1113_0085:6:2106:2276:6097_length=100  16      chr4    3000069 39      95M     *       0       0       CATTTCCAATACTATACCAAAAATCCCCCATACCTACCCACCCGCACTCCCCTACCCACCT
SRR630258.8741011_HWI-ST1113_0085:6:1208:15012:167343_length=100        0       chr4    3000448 32      95M     *       0       0       TTCGATAAAATTTTGTTAGTGTATGTAATGGTGTTAGCGTTTGGATGTTGATT
(
SRR630257.12464425_HWI-ST1113_0085:5:1308:18692:107279_length=100       16      chr4    2999991 8       95M     *       0       0       TCACCAACAAAAATTCACTCGAACAACTAAATTCTTCTTTTTTTTTTTCCATT
SRR630258.15688879_HWI-ST1113_0085:6:2106:2276:6097_length=100  16      chr4    3000069 39      95M     *       0       0       CATTTCCAATACTATACCAAAAATCCCCCATACCTACCCACCCGCACTCCCCTACCCACCT
SRR630258.8741011_HWI-ST1113_0085:6:1208:15012:167343_length=100        0       chr4    3000448 32      95M     *       0       0       TTCGATAAAATTTTGTTAGTGTATGTAATGGTGTTAGCGTTTGGATGTTGATT
SRR630257.5426488_HWI-ST1113_0085:5:1203:8069:184455_length=100 0       chr4    3001064 42      95M     *       0       0       GATTTTTAGAGTGGTTGTATAAGTTTAGAATTTTATTAATAATGGAAGAGTGTTTTTTTTT
(


cat *hapinfo.txt > Merge.hapInfo.txt
perl ../hapinfoMerge.pl Merge.hapInfo.txt

perl ~/bin/hapinfo2mhb.pl Merge.hapInfo.txt.SumUniq 0.3 > MHB.Mouse.R2-0.3.bed &
perl ~/bin/hapinfo2mhb.pl Merge.hapInfo.txt.SumUniq 0.5 > MHB.Mouse.R2-0.5.bed &
perl ~/bin/hapinfo2mhb.pl Merge.hapInfo.txt.SumUniq 0.7 > MHB.Mouse.R2-0.7.bed &
perl ~/bin/hapinfo2mhb.pl Merge.hapInfo.txt.SumUniq 0.9 > MHB.Mouse.R2-0.9.bed &


qsub SRR2011294_1.fastq.gz.job
for i in 4 5 6 7
do
qsub SRR201129$i\_1.fastq.gz.job
done


perl hapinfo2BlocAvgR2.pl mESC.sort.hapinfo.txt.SumUniq > mESC.sort.hapinfo.R2.txt
perl hapinfo2BlocAvgR2.pl Adult.sort.hapinfo.txt.SumUniq > Adult.sort.hapinfo.R2.txt

grep chr10:100004267-100004288 *sort.hapinfo.R2.txt
grep chr10:100004267-100004288 *hapinfo.txt.SumUniq
grep chr10:100004267-100004288 output.mf

grep chr10:68979553-68979781 *sort.hapinfo.R2.txt
grep chr10:68979553-68979781 *hapinfo.txt.SumUniq
grep chr10:68979553-68979781 output.mf

# R2 for each block
/home/shg047/oasis/mouse/mergeHapinfo/mESC.sort.hapinfo.R2.txt
/home/shg047/oasis/mouse/mergeHapinfo/Adult.sort.hapinfo.R2.txt
# 5me level for each block
/home/shg047/oasis/mouse/mergeHapinfo/mf/output.mf


修改hapinfo2mhb.pl的程序，输出R2,D,D-primer.
写一个程序，进行haploinfo.merge

sort Adult.hapinfo.txt > Adult.sort.hapinfo.txt
sort Adult.hapinfo.txt > Adult.sort.hapinfo.txt

data<-read.table("output.mf")
subset(data,Adult>0.6 & mESC<0.3)
 
qsub hapin2R2-mESC.job
qsub hapin2R2-adult.job


#!/bin/csh
#PBS -n hapinfo2r2
#PBS -q pdafm
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -o Adult.R2.log
#PBS -e Adult.R2.err
#PBS -V
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd /home/shg047/oasis/mouse/mergeHapinfo
perl ~/bin/hapinfo2LDR2ByBed.pl /home/shg047/oasis/mouse/RD/hapinfo/Mouse.MHB.Alice.RD90_80up_R0.5.bed Adult.hapinfo.txt > Adult.R2.txt
cd /home/shg047/oasis/mouse/mergeHapinfo
perl ~/bin/hapinfo2LDR2ByBed.pl /home/shg047/oasis/mouse/RD/hapinfo/Mouse.MHB.Alice.RD90_80up_R0.5.bed mESC.hapinfo.txt > mESC.R2.txt
cd /home/shg047/oasis/mouse/mergeHapinfo/mf
perl /home/shg047/bin/hapinfo2mf.pl /home/shg047/oasis/mouse/mergeHapinfo/mf > output.mf

sort Adult.hapinfo.txt > Adult.sort.hapinfo.txt
sort mESC.hapinfo.txt > mESC.sort.hapinfo.txt & 
perl ~/bin/hapinfo2LDR2.pl rlt chr1:10006430-10006441 < Adult.sort.hapinfo.txt


reformat.sh in=test.fq out=test.phred33.fq qin=64 qout=33

perl ~/bin/hapinfo2LDR2.pl rlt chr10:107810947-107810953  < test.sort

chr1:103067130-103067142

perl ~/bin/hapinfo2LDR2.pl rlt chr1:50329971-50335398 < Mouse.MHB.Alice.chr1.bam.hapInfo.txt
? get line-req
chr1:50329971-50335398

perl /home/shg047/bin/bam2hapInfo2PBS.pl bam2hapinfo.config  submit bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt



/home/shg047/oasis/mouse/hapinfo/Mouse.MHB.Alice.RD90_80up_R0.5.bed

/home/shg047/oasis/mouse/RD/hapinfo/



mergedBam2hapInfoV2.pl
ls SRR2011294

qsub SRR2011294

qsub SRR2011294_1.fastq.gz.job
qsub SRR2011295_1.fastq.gz.job
qsub SRR2011296_1.fastq.gz.job
qsub SRR2011297_1.fastq.gz.job


for i in `ls *cov`
do
awk '{sum +=$4;n++}END{print sum/n}' $i
done

ls -l  | awk -F : '{sum+=$5} END {print "AVG=",sum/NR}'


for i in {1..19} X Y M
do
perl ~/bin/hapinfo2mhb.pl Mouse.MHB.Alice.chr$i.bam.hapInfo.txt 0.3 >> Mouse.MHB.Alice.RD90_80up_R0.3.bed
done

hapinfo2mld.pl


perl hapinfo2mld.pl Mouse.MHB.Alice.chr$i.bam.hapInfo.txt >> Mouse.MHB.r0.25.bed

perl ~/bin/hapinfo2LDR2.pl rlt  chr14:100392580-100392941 < Mouse.MHB.Alice.chr14.bam.hapInfo.txt
samtools tview -p chr14:32546776 Mouse.MHB.Alice.chr14.bam /home/shg047/oasis/db/mm9/mm9.fa

chr14:32546776-32546898


perl /home/shg047/bin/bam2hapInfo2PBS.pl bam2hapinfo.config  submit bismark /home/shg047/oasis/db/mm9/mm9.chrom.sizes /home/shg047/oasis/db/mm9/MM9.CpG.positions.txt


for i in {1..19} X Y M
do
samtools index Mouse.MHB.Alice.chr$i.bam &
done

cd /home/shg047/oasis/db/mm9/mm9

MM9.CpG.positions.txt

/home/shg047/oasis/db/mm9/MM9.CpG.positions.txt


cat *cg.pos > MM9.CpG.positions.txt




/home/shg047/oasis/mouse/RD/MouseRD90UP80.bed

source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Mmusculus.UCSC.mm9")

chrs <- names(BSgenome.Mmusculus.UCSC.mm9)[!grepl("random", names(BSgenome.Mmusculus.UCSC.mm9))] #filter out the "upstream" chromosomes
CGs <- lapply(chrs, function(x) start(matchPattern("CG", BSgenome.Mmusculus.UCSC.mm9[[x]])))
names(CGs) <- chrs
seqlengths(genome)
genome$chr1 # same as genome[["chr1"]]


for i in {1..19} X Y M
do
awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' Mouse.MHB.Alice.chr$i.bam.RD30_80up.bed
# awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' Mouse.MHB.Alice.chr$i.bam.RD50_80up.bed
# awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' Mouse.MHB.Alice.chr$i.bam.RD70_80up.bed
# awk '{ sum += $4; n++ } END { if (n > 0) print sum; }' Mouse.MHB.Alice.chr$i.bam.RD90_80up.bed
done

for i in {1..19} X Y M
do
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' Mouse.MHB.Alice.chr$i.bam.gencov.bed 
done


awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' *bam.RD30_80up.bed
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' *bam.RD50_80up.bed
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' *bam.RD70_80up.bed
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' *bam.RD90_80up.bed




for i in {1..19} X Y M
do
cd /home/shg047/oasis/mouse/RD
awk '$4>29 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr$i.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr$i.bam.RD30.bed &
# awk '$4>49 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr$i.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr$i.bam.RD50.bed &
# awk '$4>69 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr$i.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr$i.bam.RD70.bed &
# awk '$4>89 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr$i.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr$i.bam.RD90.bed &
done

for i in {1..19} X Y M
do
cd /home/shg047/oasis/mouse/RD
awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD30.bed > Mouse.MHB.Alice.chr$i.bam.RD30_80up.bed &
#awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD50.bed > Mouse.MHB.Alice.chr$i.bam.RD50_80up.bed &
#awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD70.bed > Mouse.MHB.Alice.chr$i.bam.RD70_80up.bed &
#awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD90.bed > Mouse.MHB.Alice.chr$i.bam.RD90_80up.bed &
done



for i in {1..19} X Y M
do
cd /home/shg047/oasis/mouse/RD
awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' Mouse.MHB.Alice.chr$i.bam.RD50.bed > Mouse.MHB.Alice.chr$i.bam.RD50_80up.bed &
done





samtools view -b -q 20 Mouse.MHB.Alice.MergeBam.sort.bam chr1 > Mouse.MHB.Alice.chr1.bam &
samtools view -b -q 20 Mouse.MHB.Alice.MergeBam.sort.bam chr2 > Mouse.MHB.Alice.chr1.bam &
awk '$4>29 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr1.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr1.bam.RD10.bed &
awk '$4>29 { print $1"\t"$2"\t"$3}' Mouse.MHB.Alice.chr2.bam.gencov.bed  | bedtools merge -d 20 -i - > Mouse.MHB.Alice.chr2.bam.RD10.bed &


./configure --prefix=/home/shg047/software/R-3.3.1 '--with-cairo' \
 '--with-jpeglib' '--with-readline' '--with-tcltk' '--with-x=no'\
 '--with-blas' '--with-lapack' '--enable-R-profiling' '--with-tiff=yes'\
 '--enable-R-shlib'\
 '--enable-memory-profiling'
 
SRR1248497_1.fas
SRR630206
SRR833525_1


samtools index Mouse.MHB.Alice.MergeBam.sort.bam
for i in {1..22} X Y M
do
samtools view -b -q 20 Mouse.MHB.Alice.MergeBam.sort.bam chr$i > Mouse.MHB.Alice.chr$i.bam &
done

for i in {01,02,04,05,06,07,09,10,11,12}
do
echo $i
done

for i in {01,02,04,05,06,07,09,10,11,12}
do
cat s_*_1_Indx$i.txt.gz > Indx$i.read1.fq.gz &
cat s_*_2_Indx$i.txt.gz > Indx$i.read2.fq.gz &
done

 
zcat *1_Indx02*.gz > Indx02.s.1.fq.gz 
zcat *2_Indx02*.gz > Indx02.s.2.fq.gz 
Indx10

perl ../../fastq/bismarkPBS.pl SamConfig.txt submit

mv	s_1_Indx01.txt.gz	Indx01_s_1.txt.gz
mv	s_1_Indx02.txt.gz	Indx02_s_1.txt.gz
mv	s_1_Indx04.txt.gz	Indx04_s_1.txt.gz
mv	s_1_Indx05.txt.gz	Indx05_s_1.txt.gz
mv	s_1_Indx06.txt.gz	Indx06_s_1.txt.gz
mv	s_1_Indx07.txt.gz	Indx07_s_1.txt.gz
mv	s_1_Indx09.txt.gz	Indx09_s_1.txt.gz
mv	s_1_Indx10.txt.gz	Indx10_s_1.txt.gz
mv	s_1_Indx11.txt.gz	Indx11_s_1.txt.gz
mv	s_1_Indx12.txt.gz	Indx12_s_1.txt.gz
mv	s_2_Indx01.txt.gz	Indx01_s_2.txt.gz
mv	s_2_Indx02.txt.gz	Indx02_s_2.txt.gz
mv	s_2_Indx04.txt.gz	Indx04_s_2.txt.gz
mv	s_2_Indx05.txt.gz	Indx05_s_2.txt.gz
mv	s_2_Indx06.txt.gz	Indx06_s_2.txt.gz
mv	s_2_Indx07.txt.gz	Indx07_s_2.txt.gz
mv	s_2_Indx09.txt.gz	Indx09_s_2.txt.gz
mv	s_2_Indx10.txt.gz	Indx10_s_2.txt.gz
mv	s_2_Indx11.txt.gz	Indx11_s_2.txt.gz
mv	s_2_Indx12.txt.gz	Indx12_s_2.txt.gz


for i in {1..19} X Y M
do
cd /home/shg047/oasis/mouse/RD
awk '$4>29 { print $1"\t"$2"\t"$3}' merge.chr$i.bam.gencov.bed | bedtools merge -d 10 -i - > merge.chr$i.bamRD10.bed
awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' merge.chr$i.bamRD10.bed > merge.chr$i.bam.RD10_80up.bed
done


/home/shg047/oasis/mouse/RD
Mouse.MHB.Alice.MergeBam.bam

/oasis/tscc/scratch/zhl002/WGBS_mouse

Previous Script:
1, Download or Prepare SRA Download Configure  (http://www.ebi.ac.uk/ena/data/view/SRP028600)
2, perl fastqDownload.pl SamConfig.txt 8 submit   (fastq-dump --split-files --gzip)
3, trim_galore --phred33 --fastqc --stringency 3 -q 20 --trim1 --length 20 --gzip --clip_R1 2 --three_prime_clip_R1 2 --illumina SRR299055_1.fastq.gz --output_dir ../fastq_trim
4, bismark_genome_preparation ./     # merge all the fa to one file (mm9.fa or hg19.fa)
5, 
--------------------------------------------------------------------------------------------------------------






gawk '{ sum += $1 }; END { print sum }' file
gawk -F: '{ print $1 }' /etc/passwd

ls *_report.txt | awk -FS="_" 'BEGIN { FS = "_" } {print $1}' | sort -u | wc -l
ls *_report.txt | awk -F_ '{print $1}' | sort -u | wc -l

perl -p -i -e 's/walltime=167/walltime=48/g' *job
perl -p -i -e 's/ppn=16/ppn=4/g' *job
perl -p -i -e 's/multicore 6/multicore 2/g' *job

perl bismarkPBS.pl FastqMatchConfig.txt submit
bismark --bowtie2 --phred33-quals --fastq -L 30 -N 1 --multicore 6 /home/shg047/db/mm9 -1 ../fastq_trim/SRR833525_1.fastq.gz_val_1.fq.gz -2 ../fastq_trim/SRR833525_2.fastq.gz_val_2.fq.gz  -
qsub SRR833525_1.fastq.gz.job
qsub SRR833526_1.fastq.gz.job

/home/shg047/oasis/db/mm9
/home/shg047/oasis/mouse/aliceSRR2876131_1_val_1.fq.gz

for i in 01 02 04 05 06 07 09 10 11 12 
do 
cat s_*_1_Indx$i.txt.gz > s_1_Indx$i.txt.gz &
cat s_*_2_Indx$i.txt.gz > s_2_Indx$i.txt.gz &
done 

cat s_*_4_Indx01.txt.gz > s_4_Indx01.txt.gz &
cat s_*_5_Indx01.txt.gz > s_5_Indx01.txt.gz &
cat s_*_6_Indx01.txt.gz > s_6_Indx01.txt.gz &
cat s_*_7_Indx01.txt.gz > s_7_Indx01.txt.gz &
cat s_*_9_Indx01.txt.gz > s_9_Indx01.txt.gz &
cat s_*_10_Indx01.txt.gz > s_10_Indx01.txt.gz &
cat s_*_11_Indx01.txt.gz > s_11_Indx01.txt.gz &
cat s_*_12_Indx01.txt.gz > s_12_Indx01.txt.gz &

gzcat file1.fastq.gz file2.fastq.gz | gzip > merged.fastq.gz

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz

/media/LTS_60T/SeqStore2016/130104_SN1001
/media/LTS_60T/SeqStore2016/130104_SN1001


1, SRR2136776 require trim and alignment

/home/shg047/oasis/db/mm9

wget http://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.16.3.tar.gz
cd /home/shg047/oasis/db/mm9

# 2016-08-22
perl -p -i -e 's/walltime=7/walltime=8/g' *job
cd /home/shg047/oasis/mouse
perl fastqDownload.pl SamConfig.txt 8 submit



# 2016-08-17

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("tutorial.Rmd")'
R -e 'library("markdown");rpubsUpload("tutorial","tutorial.html")'

library("knitr")
library("rmarkdown")
pandoc("tutorial.md",format="MediaWiki")


rmarkdown::render('tutorial.Rmd')  # output: md_document


# 2016-08-16
R -e 'library("rmarkdown");library("knitr");rmarkdown::render("NormalDevconJuly.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","NormalDevconJuly.html")'
library("knitr")
library("rmarkdown")
rmarkdown::render('DeconvolutionMixture.Rmd')
pandoc("DeconvolutionMixture.md",format="MediaWiki")

1. Download geneSCF: 
2, tar zxvf 
./geneSCF -m=update -i=test/H0.list -o=test/output/ -t=sym -db=GO_MF -bg=20000 --plot=yes -org=goa_human

ls *job |cut -b 1-11 | sort -u | wc -l
# 2016-08-17
perl bam2amfBybedByPileOMethPBS.pl mhb.bed /home/shg047/oasis/N37/sortBam
perl bam2amfBybedByPileOMethPBS.pl mhb.bed /home/shg047/oasis/SALK/bam
#2016-08-12
cd /home/shg047/oasis/N37/sortBam
cd /home/shg047/oasis/SALK/bam

PileOMeth extract -q 10 -p 5 --minDepth 5 -l mhb.bed ~/oasis/db/hg19_lambda.fa /home/shg047/oasis/N37/sortBam/Indx01.sort.bam


cd /home/shg047/oasis/monod/Figure3
PileOMeth extract -q 10 -p 5 --minDepth 5 -l mhb.bed ~/oasis/db/hg19_lambda.fa
PileOMeth extract  -q 10 -p 5 --minDepth 5 /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam  -o 7-P-3.bedGraph 
#2016-08-9
data<-read.table("xx",head=T)
y<-paste(data[,5],data[,6],sep="-")
table(y)
bedtools sort -i xx.bed > xx.sort.bed
bedtools closest -fu -D "ref" -a xx.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9}' | sort -u > xx.anno.bed
# 
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=Figure-gs.pdf Figure*.pdf
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=Supp.Figure-gs.pdf Supp*Fig*.pdf

pdftk Figure*.pdf cat output Figure.pdf
pdftk Supp*Fig*.pdf cat output Supp.Figure.pdf

# 2016-07-29

bedtools sort -i HumanTissueGSI.bed > HumanTissueGSI.sort.bed
bedtools closest -fd -D ref -a HumanTissueGSI.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed


# 2016-07-21
bedtools closest -a 
http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=tlexander&hgS_otherUserSessionName=Colon.Cancer.Hg19  
http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=tlexander&hgS_otherUserSessionName=ESCA-MH450
http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=tlexander&hgS_otherUserSessionName=mm9
http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=tlexander&hgS_otherUserSessionName=MONOD

# 2016-07-20
sudo /etc/init.d/apache2 stop
sudo /etc/init.d/mysql stop
sudo /etc/init.d/proftpd stop
sudo /opt/lampp/lampp start
sudo service vsftpd stop
sudo /opt/lampp/lampp stop


cd ~/Downloads
sudo chmod 777 -R bitnami-wordpress-3.9.1-1-module-linux-x64-installer.run
./bitnami-wordpress-3.9.1-1-module-linux-x64-installer.run



install.packages('gsheet')
library("gsheet")
gsheet2tbl('https://drive.google.com/open?id=0B2TJ0NCGGrdpVlBXcWJOY2hRdFk')


screen -ls | grep pts | cut -d. -f1 | awk '{print $1}' | xargs kill

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("Deconv.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","Deconv.html")'
R -e 'library("markdown");rpubsUpload("normalDev","Deconv.html")'
R -e 'library("rmarkdown");library("knitr");rmarkdown::render("RandomForst.Rmd")'


 #!/bin/csh
 #PBS -N bedgraph2matrix
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/hapinfo/June
 R -e 'library("rmarkdown");library("knitr");rmarkdown::render("RandomForst.Rmd")'



R -e 'library("rmarkdown");library("knitr");rmarkdown::render("RandomForest.Rmd")'



---
title: "Habits"
output:
  md_document:
    variant: markdown_github
---
---
title: "Planets"
author: "Manoj Kumar"
date: "March 3, 2016"
output: v
  html_document:
    toc: true # table of content true
    depth: 3  # upto three depths of headings (specified by #, ## and ###)
    number_sections: true  ## if you want number sections at each table header
    theme: united  # many options for theme, this one is my favorite.
    highlight: tango  # specifies the syntax highlighting style
---

2016-07-11
exampl.Rmd

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("exampl.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","ColonDevconJuly.Rmd.html")'

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("ColonDevconJuly.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","ColonDevconJuly.Rmd.html")'

cmake .. -DRSTUDIO_TARGET=Server -DCMAKE_BUILD_TYPE=Release  CMAKE_INSTALL_PREFIX=./
W: GPG error: http://ppa.launchpad.net trusty InRelease: The following signatures couldn't be verified because the public key is not available: NO_PUBKEY 4F191A5A8844C542
sudo mv /var/lib/dpkg/lock
sudo mv /var/cache/apt/archives/lock 
sudo mv /var/cache/apt/archives/lock_bak
sudo rm /var/lib/apt/lists/* -vf
sudo apt-get update && sudo apt-get upgrade
sudo add-apt-repository --remove ppa:kernel-ppa/ppa

/etc/apt/sources.list
deb http://archive.ubuntu.com/ubuntu trusty main
lsb_release -sc
sudo apt-get update
sudo apt-get install gedbi-core

R -e 'library("rmarkdown");library("knitr");rmarkdown::render("NormalDevconJuly.Rmd")'
R -e 'library("markdown");rpubsUpload("normalDev","NormalDevconJuly.html")'

R CMD BATCH RunchampDMR.R &
track type=bigWig name=proteinA smoothingWindow=4 color=123,100,50 autoScale=on viewLimits=1:200 visibility=full windowingFunction=maximum bigDataUrl=https://projects/files/file.bw
48502
chr1:23730471-23730518
#!/usr/bin/sh
for i in `ls *.bedGraph`
do 
(awk 'NR!=1' $i| sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4}') > $i.dehsort 
bedGraphToBigWig $i.dehsort ~/oasis/db/hg19.chrom.sizes $i.dehsort.bw 
bigWigAverageOverBed $i.dehsort.bw  /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed  $i.dehsort.bw.mf
rm *.dehsort
done
head -n 20 /home/shg047/oasis/monod/BAM/MF_PileOMeth/CRC-P-002_CpG.bedGraph
head -n 20 /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/6-P-2.sorted.clipped_CpG.bedGraph
 samtools tview /home/shg047/oasis/monod/BAM/rename/CRC-P-024.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:23730471-23730518
cp /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/mf2matrix.pl  ./
cp /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/bedgraph2matrix.pl ./
cp /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/bw2mf.pl ./
 
 
 
 samtools tview RRBS-7P23.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr10:100027957
samtools tview PC-P-5 /home/shg047/oasis/db/hg19/hg19.fa -p chr1:17092527-17092562
samtools tview PC-P-5 /home/shg047/oasis/db/hg19/hg19.fa -p 
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:10496-10498

samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:935977-936090
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066


 
 #!/bin/csh
 #PBS -N bedgraph2matrix
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

perl bedgraph2matrix.pl > bedgraphMatrix.txt
 
 

 
 
 

 

 #!/bin/csh
 #PBS -N deconv
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/hapinfo/June
 #Rscript --vanilla colon.updateGSI.R
 #Rscript --vanilla lung.updateGSI.R
 Rscript --vanilla normal.updateGSI.R

 
 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

(tail -n +2 6-T-4.sorted.clipped_CpG.bedGraph| sort -k 1,1 -k2,2n) > 6-T-4.sorted.clipped_CpG.bedGraph


(awk 'NR!=1' 6-T-4.sorted.clipped_CpG.bedGraph| sort -k 1,1 -k2,2n | awk '{print $1,$2,$3,$4}') > 6-T-4.sorted.clipped_CpG.sort.bedGraph
bedGraphToBigWig 6-T-4.sorted.clipped_CpG.sort.bedGraph ~/oasis/db/hg19.chrom.sizes 6-T-4.sorted.clipped_CpG.sort.bedGraph.bw

 
 

# copy simulated CCTmix and LCTmix to hapinfo fold, run hapinfo2mhl perl script
cp
cp 

cd /home/shg047/oasis/monod/hapinfo/hapinfo
ls *hapInfo.txt > Hapinfo_File_list

qsub hapinfo2mhl.pbs




perl CCTmixturePBS.pl 6-P-Merge.hapinfo.txt WB.hapInfo.txt
perl LCTmixturePBS.pl 7-P-Merge.hapinfo.txt WB.hapInfo.txt
  
echo "perl /home/shg047/bin/hapinfo2mhl.pl Hapinfo_File_list > MHL.output.txt" |qsub -N hapinfo2MHL -q pdafm -l walltime=8:00:00
echo "perl /home/shg047/bin/hapinfo2mhl.pl Hapinfo_File_list > MHL.output.txt" |qsub -N h2m

hapInfo.txt

 #!/bin/csh
 #PBS -N bam2MHB
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group 
 perl /home/shg047/bin/hapinfo2mhl.pl Hapinfo_File_list > MHL.mixture.simulation.txt


#!/usr/bin/env perl
use strict;
my @f=qw /0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/;
foreach my $f(@f){
my $cmd1="perl ../randomSampleFromHaploInfo.pl int(24850020*f) WB.hapinfo.merge > WB";
my $cmd2="perl ../randomSampleFromHaploInfo.pl int(1386115*(1-f)) 6-T-Merge.hapinfo.txt > 6T";
my $cmd3="cat WB 6T > Wcmix.f.hapinfo";
print "$cmd1\n$cmd2\n$cmd3\n";
}

/home/shg047/oasis/monod/deconvolution/simulation
# 2016-06-21

 #!/bin/csh
 #PBS -N bam2MHB
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 Rscript --vanilla /home/shg047/oasis/monod/hapinfo/colon.updateGSI.R
 Rscript --vanilla /home/shg047/oasis/monod/hapinfo/lung.updateGSI.R
 Rscript --vanilla /home/shg047/oasis/monod/hapinfo/lungcancer.updateGSI.withoutLCT.R
 Rscript --vanilla /home/shg047/oasis/monod/hapinfo/coloncancer.updateGSI.withoutCCT.R

echo "Rscript --vanilla /home/shg047/oasis/monod/hapinfo/colon.updateGSI.R" | qsub -q pdafm
echo "Rscript --vanilla /home/shg047/oasis/monod/hapinfo/lung.updateGSI.R" | qsub -q pdafm
echo "Rscript --vanilla /home/shg047/oasis/monod/hapinfo/lungcancer.updateGSI.withoutLCT.R" | qsub -q pdafm
echo "Rscript --vanilla /home/shg047/oasis/monod/hapinfo/coloncancer.updateGSI.withoutCCT.R" | qsub -q pdafm
perl -p -i -e 's/\>0.01/\>0.05/' colon.updateGSI.R
perl -p -i -e 's/\>0.01/\>0.05/' lung.updateGSI.R
perl -p -i -e 's/\>0.01/\>0.05/' lungcancer.updateGSI.withoutLCT.R
perl -p -i -e 's/\>0.01/\>0.05/' coloncancer.updateGSI.withoutCCT.R
perl format.pl lsMHL.cor.txt | bedtools sort > tsMHL.cor.sort.bed
perl format.pl lsMHL.cor.txt | bedtools sort >lsMHL.cor.sort.bed
bedtools closest -fu -D "ref" -a lsMHL.cor.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$11}' | sort -u > lsMHL.cor.sort.anno.bed
perl format.pl php.bed |bedtools sort > php.sort.bed
bedtools closest -fu -D "ref" -a php.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$9}' | sort -u > php.bed.anno.bed
bedtools intersect -a php.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$9}' | sort -u > php.bed.anno.bed
perl format.pl | bedtools sort > tsMHL.cor.sort.bed


bedtools closest -fu -D "ref" -a WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9}' | sort -u > WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.sort.gene.anno.bed
bedtools closest -fu -D "ref" -a tsMHL.cor.sort.bed -b ~/oasis/db/hg19.refGene.sorted.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9}' | sort -u > tsMHL.cor.sort.anno.bed

# 2016-06-13

system("wget https://cran.r-project.org/src/contrib/pbkrtest_0.4-6.tar.gz")
system("wget https://cran.r-project.org/src/contrib/quantreg_5.26.tar.gz")
system("wget https://cran.r-project.org/src/contrib/lme4_1.1-12.tar.gz")
system("wget https://cran.r-project.org/src/contrib/minqa_1.2.4.tar.gz")
system("wget https://cran.r-project.org/src/contrib/nloptr_1.0.4.tar.gz")
system("wget https://cran.r-project.org/src/contrib/RcppEigen_0.3.2.8.1.tar.gz")
system("wget https://cran.r-project.org/src/contrib/SparseM_1.7.tar.gz")
system("wget https://cran.r-project.org/src/contrib/MatrixModels_0.4-1.tar.gz")
system("wget https://cran.r-project.org/src/contrib/ggplot2movies_0.0.1.tar.gz")
install.packages("minqa_1.2.4.tar.gz")
install.packages("nloptr_1.0.4.tar.gz")
install.packages("RcppEigen_0.3.2.8.1.tar.gz")
install.packages("lme4_1.1-12.tar.gz")
install.packages("pbkrtest_0.4-6.tar.gz")
install.packages("SparseM_1.7.tar.gz")
install.packages("MatrixModels_0.4-1.tar.gz")
install.packages("quantreg_5.26.tar.gz")
install.packages("car_2.1-2.tar.gz")
install.packages("ggplot2movies_0.0.1.tar.gz")

mix="simulation.mixture.inp.txt"
pure="simulation.signature.inp.txt"
class="simulation.lab.inp.txt"
Rscript fix_matrices.R $pure $mix $class
#R CMD ~/R/x86_64-pc-linux-gnu-library/3.2/Rserve/libs/Rserve --no-save
R CMD /home/shicheng/R/x86_64-pc-linux-gnu-library/3.2/Rserve/libs --no-save
java -Xmx3g -Xms3g -jar CIBERSORT.jar -M $mix.mod -P $pure.mod -c $class > outputs


# 2016-06-13
cd /home/shg047/deconvolution/cibersort
java -Xmx3g -Xms3g -jar /home/shg047/software/CIBERSORT/CIBERSORT.jar -M 

R CMD /home/shg047/software/R-3.3.0/library/Rserve/R/Rserve --no-save
R CMD /home/shg047/software/R-3.3.0/library/Rserve/libs/Rserve --no-save
R CMD ~/R/x86_64-pc-linux-gnu-library/3.2/Rserve/libs/Rserve --no-save

java -Xmx3g -Xms3g -jar /home/shg047/software/CIBERSORT/CIBERSORT.jar -M 

cp signatures.lung.data.txt.mod ~
cp test.lung.data.txt.mod ~


java -Xmx3g -Xms3g -jar /home/ddiep/softwares/CIBERSORT/CIBERSORT.jar -M colon.mixture.inp.txt -B colon.signature.inp.txt.mod

# 2016-06-13
the difference between nucleo and whole-cell 

perturbation screen 
genotype to gene expression profile
how to select postive cells (50% percent)
permuation to remove false postive gene set and check the remain set in the comparison





## Use IGV in Windows connecting with Linux server. 

In windows, you don't have terminal, you can use putty to connect Linux (installing IGV) server.
Be careful, You can not use Winscp to use Putty, since you can not set X11 forwarding for putty in Winscp. 

# before connect to linux sever with putty, please install and open Xming (Windows)

in the windows: install winscp and Xming
X11 forwarding in putty: https://wiki.utdallas.edu/wiki/display/FAQ/X11+Forwarding+using+Xming+and+PuTTY
in the linxu server install IGV for linux
Download (Binary Distribution): http://www.broadinstitute.org/igv/download
open Xming in windowns and Go to server IGV fold and run: java -Xmx750m -jar igv.jar


## Use IGV in Linux or Mac connecting with Linux server. 

ssh -X shg047@genome-miner.ucsd.edu


# 2016-06-01

/home/shg047/oasis/WGBSDB/PBMC/WGBS.select.1501.txt



http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer8/tracks_hg19/

# PBMC and Blood cells
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_BCell/tracks_hg19/Human_BCell.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_CD133HSC/tracks_hg19/Human_CD133HSC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_HSPC/tracks_hg19/Human_HSPC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_Neut/tracks_hg19/Human_Neut.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_PBMC/tracks_hg19/Human_PBMC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-100yr/tracks_hg19/Human_CD4T-100yr.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-Newborn/tracks_hg19/Human_CD4T-Newborn.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Li-PBMC-2010/Human_PBMC/tracks_hg19/Human_PBMC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-DNMT3BMut-2012/Human_BCell-Healthy/tracks_hg19/Human_BCell-Healthy.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_Macrophage/tracks_hg19/Human_Macrophage.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_Tcell/tracks_hg19/Human_Tcell.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_NK/tracks_hg19/Human_NK.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_HSC/tracks_hg19/Human_HSC.meth.bw & 
# Lung
wget http://smithlab.usc.edu/methbase/data/Xie-Human-2013/Human_IMR90/tracks_hg19/Human_IMR90.meth.bw
wget https://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Lung/UCSD.Lung.Bisulfite-Seq.STL002.wig.gz
wget https://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/IMR90_Cell_Line/UCSD.IMR90.Bisulfite-Seq.combined.wig.gz

gunzip *.gz

chr10:60515-71797

# Colon
wget http://smithlab.usc.edu/methbase/data/Berman-Human-2012/Human_ColonCancer/tracks_hg19/Human_ColonCancer.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Ziller-Human-2013/Human_Colon_Tumor_Primary/tracks_hg19/Human_Colon_Tumor_Primary.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hansen-Human-2011/Human_ColonCancer/tracks_hg19/Human_ColonCancer.meth.bw &



# Total 45 samples 
https://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Adipose_Tissue/UCSD.Adipose_Tissue.Bisulfite-Seq.STL003.wig.gz

# bedgraph to bigwig
cd /home/shg047/oasis/Estellar2016/MF
for i in `ls *bedGraph`
do
(head -n 1 $i && tail -n +2 $i | sort -k1,1 -k2,2n | awk '{print $1,$2,$3,$4}' OFS="\t" )  > $i.sort
bedGraphToBigWig $i.sort ~/oasis/db/hg19/hg19.chrom.sizes $i.bw 
done


# bigWigToBedGraph
for i in `ls *bw`
do
bigWigToBedGraph $i $i.bg &
done

# bigWigAverageOverBed
for i in `ls *bw`
do
bigWigAverageOverBed $i /home/shg047/oasis/WGBSDB/PBMC/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.renew.bed $i.aob &
done

chr18:74692045–74692213
chr10:79936919–79937042




# 2016-05-29

[shg047@tscc-login1 MF_PileOMeth]$ wc -l CTR106.bedGraphV4
20653494 CTR106.bedGraphV4
[shg047@tscc-login1 MF_PileOMeth]$ wc -l CTR107.bedGraphV4
22745883 CTR107.bedGraphV4



cd /home/shg047/oasis/DennisLo2015/MF_PileOMeth

cp /home/shg047/oasis/monod/wig/bedGraph2V4.pl /home/shg047/oasis/DennisLo2015/MF_PileOMeth

CTR101_trimmed.fq.gz_bismark_bt2.sort_CpG.bedGraph

track type="bedGraph" description="/home/shg047/oasis/DennisLo2015/MF_PileOMeth/CTR108_trimmed.fq.gz_bismark_bt2.sort CpG merged methylation levels"
chr10   61056   61058   0       0       1
chr10   61332   61334   100     1       0
chr10   61652   61654   100     2       0
chr10   63584   63586   75      3       1
chr10   64773   64775   100     9       0

vim bedGraph2V4.pl





# check what happened to the MF and MHL result

chr10:100027918-100027944       TTTTTTT 2       100027918,100027922,100027925,100027927,100027930,100027938,100027944
chr10:100027957-100027992       C       1       100027957
chr10:100027957-100027992       T       1       100027957
chr10:100121047-100121066       TTT     1       100121047,100121049,100121066
chr10:100121047-100121066       CCC     3       100121047,100121049,100121066
chr10:100121047-100121066       TT      1       100121047,100121049
chr10:100151033-100151068       T       2       100151068
chr10:100206398-100206510       TTTTTT  1       100206398,100206407,100206455,100206475,100206479,100206487
chr10:100227616-100227677       TTTT    3       100227616,100227620,100227625,100227651
chr10:100227616-100227677       TTTTTT  1       100227616,100227620,100227625,100227651,100227666,100227677
chr10:100227698-100227747       TTTTTTTT        1       100227698,100227707,100227710,100227714,100227719,100227724,100227727,100227740
chr10:100227781-100227823       TT      1       100227814,100227823
chr10:100227781-100227823       TTT     1       100227781,100227787,100227814
chr10:100992364-100992404       CT      1       100992364,100992404
chr10:100992364-100992404       C       3       100992364
chr10:100995252-100995487       T       1       100995487
chr10:100995942-100996028       T       3       100996028
chr10:101089155-101089163       TTTT    1       101089155,101089157,101089160,101089163
chr10:101089204-101089305       TTTT    2       101089204,101089209,101089243,101089250
chr10:101089204-101089305       TTTTTTTCTTT     1       101089204,101089209,101089243,101089250,101089254,101089257,101089260,101089279,101089284,101089303,101089305
chr10:101089204-101089305       TTTTCTTT        2       101089250,101089254,101089257,101089260,101089279,101089284,101089303,101089305
chr10:101089382-101089424       T       3       101089382
chr10:101089443-101089461       TTTT    2       101089443,101089454,101089456,101089461
chr10:101089481-101089508       T       2       101089481
chr10:101089481-101089508       TTTTTT  1       101089481,101089485,101089491,101089494,101089506,101089508
chr10:101089526-101089565       TTTT    1       101089532,101089537,101089554,101089565
chr10:101089526-101089565       T       1       101089526
chr10:101089526-101089565       TTTTT   2       101089526,101089532,101089537,101089554,101089565
chr10:101089569-101089596       TTT     3       101089569,101089591,101089596
chr10:101089668-101089771       TTTTTT  1       101089746,101089759,101089761,101089765,101089769,101089771
chr10:101089816-101089850       TTTTT   2       101089816,101089820,101089824,101089829,101089850
chr10:101089870-101089906       TTTT    1       101089870,101089886,101089895,101089906
chr10:101089870-101089906       T       1       101089906
chr10:101089870-101089906       TTT     1       101089870,101089886,101089895
chr10:101089928-101089935       TTT     2       101089928,101089932,101089935
chr10:101089949-101089965       TTTTT   2       101089949,101089951,101089954,101089956,101089965
chr10:101190105-101190169       TTTT    1       101190105,101190120,101190141,101190169
chr10:101190105-101190169       TTTTT   2       101190105,101190120,101190141,101190157,101190169
chr10:101190297-101190505       TTT     1       101190297,101190310,101190324
chr10:101190297-101190505       TTTTTTTTTTTTTTTT        1       101190297,101190310,101190324,101190420,101190423,101190429,101190435,101190437,101190451,101190456,101190458,101190475,101190477,101190489,101190502,101190505
chr10:101190297-101190505       TTTTTTTTTTTTT   2       101190420,101190423,101190429,101190435,101190437,101190451,101190456,101190458,101190475,101190477,101190489,101190502,101190505
chr10:101190297-101190505       TTT     1       101190384,101190389,101190420
chr10:101190297-101190505       TTTT    1       101190360,101190384,101190389,101190420
chr10:101190297-101190505       TTTTTTTTTTT     1       101190420,101190423,101190429,101190435,101190437,101190451,101190456,101190458,101190475,101190477,101190489
chr10:101190527-101190864       TTTTTTTTTTTTTTT 1       101190595,101190604,101190610,101190619,101190635,101190646,101190674,101190691,101190695,101190699,101190702,101190727,101190734,101190751,101190763
chr10:101190527-101190864       TTTTTTTTTT      1       101190527,101190558,101190565,101190595,101190604,101190610,101190619,101190635,101190646,101190674
chr10:101190527-101190864       TTTTT   1       101190776,101190798,101190801,101190815,101190827
chr10:101190527-101190864       TT      1       101190841,101190864
chr10:101190527-101190864       TTT     4       101190527,101190558,101190565


RRBS-7P23.hapInfo.txt
ls RRBS-7P23*
less RRBS-7P23.sorted.clipped_CpG.bedGraph

grep 540439 RRBS-7P23.hapInfo.txt | grep chr1
grep 101190298 RRBS-7P23.sorted.clipped_CpG.bedGraph
grep 101190120 RRBS-7P23.sorted.clipped_CpG.bedGraph
grep 101190120 RRBS-7P23.sorted.clipped_CpG.bedGraph



less -S RRBS-7P23.hapInfo.txt

less 

samtools tview RRBS-7P23.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr10:100027957
samtools tview PC-P-5 /home/shg047/oasis/db/hg19/hg19.fa -p chr1:17092527-17092562
samtools tview PC-P-5 /home/shg047/oasis/db/hg19/hg19.fa -p 
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:935977-936090
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066
samtools tview ../bam/PC-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:876056-876066





samtools tview RRBS-7P23.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:850885


chr1    850885

chr10:100027957-100027992

/home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/6-P-9.sorted.clipped_CpG.bedGraph
wc -l 6-P-9_CpG.bedGraph

cd oasis/monod/bam/RRBS1/bam/

chr1:10524-10526
samtools tview 6-P-9.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:10524
samtools tview 6-P-9.sorted.clipped.bam /home/shg047/oasis/db/hg19/hg19.fa -p chr1:724245

chr1:724245


wc -l ../../../wig/6-P-9_CpG.bedGraph
1178863

head ../../../bedGraph/6-P-9_CpG.bedGraph

head 6-P-9_CpG.bedGraph
1394266 6-P-9.sorted.clipped_CpG.bedGraph





http://smithlab.usc.edu/methbase/data/Hon-Human-2012/Human_HCC1954/tracks_hg19/Human_HCC1954.meth.bw
# 2016-05-27
Install BioPerl Without Root Privileges in Ubuntu/Linxu
https://www.biostars.org/p/193668/

< How To Install BioPerl Without Root Privileges in Ubuntu/Linxu>

# Install certain nessary library
perl -MCPAN -Mlocal::lib -e 'CPAN::install(local::lib)'

# Download latest bioperl
git clone https://github.com/bioperl/bioperl-live.git
cd bioperl-live
perl Build.PL
./Build test 

# biuld test wrong?? then test them one by on, such as 
./Build test --test_files t/LocalDB/Taxonomy/sqlite.t --verbose
./Build test --test_files t/Root/RootIO.t --verbose

# install the failed library according bulid test result
perl -MCPAN -Mlocal::lib -e 'CPAN::install(DBI)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(DBD::SQLite)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(LWP::UserAgent)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(LWP::UserAgent)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::Phylo)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::DB::Sam)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Graph::Directed)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Twig)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::Ext::Align)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Parser::PerlSAX)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Simple)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(HTML::TableExtract)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(IO::Scalar)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::SAX)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::SeqIO::staden::read)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Parser::PerlSAX)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::SeqIO::staden::read)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::LibXML)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Convert::Binary::C)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::SAX::Writer)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(XML::Writer)'
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::SeqIO::staden::read)'

perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::DB::Sam)'

Bio::DB::Sam
perl -MCPAN -Mlocal::lib -e 'CPAN::install(Bio::DB::Sam)'

# test again
./Build test 
## Test bioperl installation (you should get a version number)
perl -MBio::Root::Version -le 'print $Bio::Root::Version::VERSION'



wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2
tar xjf samtools-0.1.18.tar.bz2 && cd samtools-0.1.18
make CFLAGS=-fPIC
export SAMTOOLS=`pwd`
cpanm Bio::DB::Sam




# 2016-05-26

cd /home/shg047/oasis/SALK/bam
# only call SNP
samtools mpileup --skip-indels -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f <reference genome.fa> -o <output.bcf> <list of input bam files>
bcftools index  <output.bcf> <indexed.bcf>
bcftools call --skip-variants indels --multiallelic-caller --variants-only  -O v <output.bcf> -o <output.vcf>
# only call Indel
samtools mpileup -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f <reference genome.fa> -o <output.bcf> <list of input bam files>
bcftools index  <output.bcf> <indexed.bcf>
bcftools call --skip-variants snps --multiallelic-caller --variants-only  -O v <output.bcf> -o <output.vcf>
# call indel and SNP
samtools mpileup -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f <reference genome.fa> -o <output.bcf> <list of input bam files>
bcftools index  <output.bcf> <indexed.bcf>
bcftools call --multiallelic-caller --variants-only  -O v <output.bcf> -o <output.vcf>

samtools mpileup -d 250 -m 1 -E --BCF --output-tags DP,AD,ADF,SP -f /home/shg047/oasis/db/hg19/hg19.fa -o STL001RV-01.chr6.sorted.clipped.bcf  STL001RV-01.chr6.sorted.clipped.bam
bcftools index  STL001RV-01.chr6.sorted.clipped.bcf  
bcftools call --skip-variants indels --multiallelic-caller --variants-only  -O v <output.bcf> -o <output.vcf>

bcftools call --multiallelic-caller --variants-only  -O v STL001RV-01.chr6.sorted.clipped.bcf  -o STL001RV-01.chr6.sorted.clipped.vcf
bcftools call --skip-variants snps --multiallelic-caller --variants-only  -O v STL001RV-01.chr6.sorted.clipped.bcf  -o STL001RV-01.chr6.sorted.clipped.vcf

cd  /oasis/tscc/scratch/ddiep/Working/160428_RapidRun
You must learn the pipeline for the methylation calling

rm * chrLambdaNEB.methylFreq

cd /home/shg047/oasis/monod/hapinfo
R CMD BATCH trim.R

/home/shg047/oasis/monod/HMHPlasma/Excl/excl.txt
742746b3f183
xse251<-unique(xse)
match(xse251,xse)

Writing reviews of academic papers: https://github.com/jtleek/reviews
A guide to reading scientific papers: GitHub - jtleek/readingpapers: 
Good Habit for Bioinformatics Analyst or Scientist: https://www.biostars.org/p/190366/
Lab notebook for bioinformatics/data analysis work: https://www.biostars.org/p/191984/
The Most Common Stupid Mistakes In Bioinformatics: https://www.biostars.org/p/7126/#191667

cp ../Bladder.hapInfo.txt ./ 
cp ../Brain.hapInfo.txt ./ 
cp ../Colon.hapInfo.txt ./ 
cp ../Esophagus.hapInfo.txt ./ 
cp ../Intestine.hapInfo.txt ./ 
cp ../Kidney.hapInfo.txt ./ 
cp ../Liver.hapInfo.txt ./ 
cp ../Lung.hapInfo.txt ./ 
cp ../Pancreas.hapInfo.txt ./ 
cp ../Stomach.hapInfo.txt ./ 
cp ../NC-P.hapInfo.txt ./ 
cp ../WB.hapInfo.txt ./ 


#!/usr/bin/perl
use CWD;
my $dir=getcwd;
my @file=glob("NC-P*hapInfo.txt");
foreach my $file(@file){
my ($sam,undef)=split/\./,$file;
open OUT,">$file.job";
print OUT "#!/bin/csh\n";
print OUT "#PBS -n $file\n";
print OUT "#PBS -q hotel\n";
print OUT "#PBS -l nodes=1:ppn=1\n";
print OUT "#PBS -l walltime=1:00:00\n";
print OUT "#PBS -V\n";
print OUT "#PBS -M shicheng.guo\@gmail.com \n";
print OUT "#PBS -m abe\n";
print OUT "#PBS -A k4zhang-group\n";
print OUT "cd $dir\n";
print OUT "perl HMHInPlasmaTest.pl $file excl.txt $sam.HMH\n";
}

system("wget https://cran.r-project.org/src/contrib/gdata_2.17.0.tar.gz")
system("wget https://cran.r-project.org/src/contrib/caTools_1.17.1.tar.gz")

Step 1. Download Roadmap Data
cd /home/shg047/oasis/Roadmap
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/EG.mnemonics.name.xls
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/FractionalMethylation.tar.gz
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/FractionalMethylation.tar.gz.md5sum
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/header
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/ReadCoverage.tar.gz
wget ftp://ftp.bcgsc.ca/public/mbilenky/112epigenomes/5mC/SBS_Removed_E027_E064_Fixed_E012/ReadCoverage.tar.gz.md5sum

Step 1. Repalce -1 to space, so that perl could calculate the mean and SD based on that table
for i in {1..22} "X" "Y" "M"
do
perl -p -i -e "s/-1/ /g" chr$i.fm &
done

Step 2. Calculate the Regionss mean methylation level and Summary SD

perl MethSearchRoadmap.pl target.txt ESCA-space.fm.rlt.txt 23 2 &
perl MethSearchRoadMapByLoc.pl target.txt Target-Full-Tissue-space.fm.rlt.txt 2 &

wget http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/DMRs/WGBS_DMRs_v2.tsv.gz

for i in {1..5}
do
echo "perl /home/shg047/oasis/monod/HMHPlasma/HMHInPlasmaTest.pl 6-P-1.hapInfo.txt 6-T-1.hapInfo.txt excl.txt 6-P-1.HMH" | qsub -q hotel
echo "perl /home/shg047/oasis/monod/HMHPlasma/HMHInPlasmaTest.pl 7-P-1.hapInfo.txt 7-T-1.hapInfo.txt excl.txt 7-P-1.HMH" | qsub -q hotel
echo "perl /home/shg047/oasis/monod/HMHPlasma/HMHInPlasmaTest.pl PC-P-1.hapInfo.txt PC-T-1.hapInfo.txt excl.txt PC-P-1.HMH" | qsub -q hotel
done

 #!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test
 #PBS -e test
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 
perl /home/shg047/oasis/monod/HMHPlasma/HMHInPlasmaTest.pl 6-P-1.hapInfo.txt 6-T-1.hapInfo.txt excl.txt 6-P-1.HMH
perl ./HMHInPlasmaTest.pl 7-P-1.hapInfo.txt 7-T-1.hapInfo.txt excl.txt 7-P-1.HMH
perl ./HMHInPlasmaTest.pl PC-P-1.hapInfo.txt PC-T-1.hapInfo.txt excl.txt PC-P-1.HMH
perl ./HMHInPlasmaTest.pl 6-P-2.hapInfo.txt 6-T-2.hapInfo.txt excl.txt 6-P-2.HMH
perl ./HMHInPlasmaTest.pl 7-P-2.hapInfo.txt 7-T-2.hapInfo.txt excl.txt 7-P-2.HMH
perl ./HMHInPlasmaTest.pl PC-P-2.hapInfo.txt PC-T-2.hapInfo.txt excl.txt PC-P-2.HMH
perl ./HMHInPlasmaTest.pl 6-P-3.hapInfo.txt 6-T-3.hapInfo.txt excl.txt 6-P-3.HMH
perl ./HMHInPlasmaTest.pl 7-P-3.hapInfo.txt 7-T-3.hapInfo.txt excl.txt 7-P-3.HMH
perl ./HMHInPlasmaTest.pl PC-P-3.hapInfo.txt PC-T-3.hapInfo.txt excl.txt PC-P-3.HMH
perl ./HMHInPlasmaTest.pl 6-P-4.hapInfo.txt 6-T-4.hapInfo.txt excl.txt 6-P-4.HMH
perl ./HMHInPlasmaTest.pl 7-P-4.hapInfo.txt 7-T-4.hapInfo.txt excl.txt 7-P-4.HMH
perl ./HMHInPlasmaTest.pl PC-P-4.hapInfo.txt PC-T-4.hapInfo.txt excl.txt PC-P-4.HMH
perl ./HMHInPlasmaTest.pl 6-P-5.hapInfo.txt 6-T-5.hapInfo.txt excl.txt 6-P-5.HMH
perl ./HMHInPlasmaTest.pl 7-P-5.hapInfo.txt 7-T-5.hapInfo.txt excl.txt 7-P-5.HMH
perl ./HMHInPlasmaTest.pl PC-P-5.hapInfo.txt PC-T-5.hapInfo.txt excl.txt PC-P-5.HMH


#!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test
 #PBS -e test
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/Batch3Plasma/hapinfo
 perl ~/bin/hapinfo2mf.pl ./ > batch3AMF.txt
 perl ~/bin/hapinfo2mhl.pl  Hapinfo_File_list > batch3MHL.txt

qsub NC-P-43.*job
qsub CRC-P-3.*job
qsub CRC-P-10.*job
qsub NC-P-47.*job



cd  /oasis/tscc/scratch/ddiep/Plasma_10ngRRBS/BAMfiles
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl > /home/shg047/oasis/monod/Batch3Plasma/saminfo.txt
cd /home/shg047/oasis/monod/Batch3Plasma/
cp *bam *bai /home/shg047/oasis/monod/Batch3Plasma/bam

perl ~/bin/bam2hapInfo2PBS.pl ../saminfo.txt submit bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt


7-P-1.hapInfo.txt
a='80981116-80981175'

grep $a 6-P-1.hapInfo.txt
grep $a 6-T-1.hapInfo.txt
grep $a NC-P.hapInfo.txt



cp ../hapinfo/6-P-1.hapInfo.txt ./
cp ../hapinfo/6-P-2.hapInfo.txt ./
cp ../hapinfo/6-P-3.hapInfo.txt ./
cp ../hapinfo/6-P-4.hapInfo.txt ./
cp ../hapinfo/6-P-5.hapInfo.txt ./
cp ../hapinfo/7-P-1.hapInfo.txt ./
cp ../hapinfo/7-P-2.hapInfo.txt ./
cp ../hapinfo/7-P-3.hapInfo.txt ./
cp ../hapinfo/7-P-4.hapInfo.txt ./
cp ../hapinfo/7-P-5.hapInfo.txt ./
cp ../hapinfo/PC-P-1.hapInfo.txt ./
cp ../hapinfo/PC-P-2.hapInfo.txt ./
cp ../hapinfo/PC-P-3.hapInfo.txt ./
cp ../hapinfo/PC-P-4.hapInfo.txt ./
cp ../hapinfo/PC-P-5.hapInfo.txt ./

perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-2.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-3.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1
perl report_cancerHM.pbs 6-T-1.hapInfo.txt 6-T-1.hapInfo.txt NC-P.hapInfo.txt 6-T-1


report.CanHMH.pl

/home/shg047/bak/plink/china

cd /oasis/tscc/scratch/shg047/monod/hapinfo
head -n 2000 6-P-1.hapInfo.txt > /oasis/tscc/scratch/shg047/monod/test/6-P-1.hapInfo.txt
head -n 2000 6-P-1.hapInfo.txt > /oasis/tscc/scratch/shg047/monod/test/6-T-1.hapInfo.txt

cat NC-P-*.hapInfo.txt >> /oasis/tscc/scratch/shg047/monod/test/NC-P.hapInfo.txt
cd /oasis/tscc/scratch/shg047/monod/test/


wig2bed  --zero-indexed  < wgEncodeBroadHistoneGm12878H3k4me1StdSig.wig  > wgEncodeBroadHistoneGm12878H3k4me1StdSig.bed
bedtools cluster -d 50 -i wgEncodeBroadHistoneGm12878H3k4me1StdSig.wig > wgEncodeBroadHistoneGm12878H3k4me1StdSig.bed.bed &
./PileOMeth extract -r chr1:723205-725116 --minDepth 1 /home/shg047/oasis/db/hg19/hg19.fa /home/shg047/oasis/monod/bam/RRBS1/bam/PC-T-7.sorted.clipped.bam -o /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/PC-T-7.sorted.clipped
./PileOMethG extract -r chr1:723205-725116 --minDepth 1 /home/shg047/oasis/db/hg19/hg19.fa /home/shg047/oasis/monod/bam/RRBS1/bam/PC-T-7.sorted.clipped.bam -o /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/PC-T-7.sorted.clipped

SIEMENS
awk '{if ($0 ~ /^>/) { print $0; } else { printf("%s%c%s\n", substr($0, 1, 9), "X", substr($0, 11, length($0) - 10))}}' in.fa > out.fa
  
/home/shg047/oasis/db/mm9/chr10.fa
samtools tview Indx01.merged.bam.sorted.bam /home/shg047/oasis/db/mm9/chr10.fa 

chr10:100222302-100223142


samtools bedcov  /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  /home/shg047/oasis/monod/bam/RRBS2/RRBS2/RRBS-6P16.sorted.clipped.bam 


grep chr10:10003965-10004290 B6ES.hapInfo_fixed.txt

grep chr10:100222148-100223142 2A4F1_miPS.hapInfo_fixed.txt

samtools tview 

chr10:100222148-100223142

#!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test
 #PBS -e test
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /oasis/tscc/scratch/shg047/monod/test
 perl report_cancerHM.pl ../hapinfo/6-P-1.hapInfo.txt ../hapinfo/6-T-1.hapInfo.txt ../hapinfo/NC-P-1.hapInfo.txt
 
 
 cd /home/shg047/bak/plink/phase  
 phase -f0 -F0.1 -p0.7 -c test.inp test.out


 cd /oasis/tscc/scratch/zhl002/hapinfo
 perl ~/bin/hapinfo2mhl.pl Hapinfo_File_List.txt > ~/MHL.output.12april.txt

 cd /oasis/tscc/scratch/zhl002/hapinfo

 grep chr10:100017509-100017813 B6ES.hapInfo_fixed.txt
 

2016-04-11

grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/7-T-2*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/7-T-3*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/7-P-9*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/7-T-1*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/DennisLo2015/hapinfo/CTR118*hapInfo.txt
grep chr19:58220626-58220668 /home/shg047/oasis/monod/hapinfo/STL003SG-01*hapInfo.txt

7-T-1
torrow:

1) transfer 6 column file to 4 column file
2) try to upload to UCSC
3, try to find what happend to the plasma data?? any thing wrong? and check the row data



 #!/usr/bin/perl
 foreach my $file (glob("*bw")){
 my ($name,undef)=split /\_/,$file;
 print "track type=bigWig name=\"$name\" description=\"$file\" bigDataUrl\=ftp\:\/\/ucsd002\:ucsd002\@132.239.189.199\/wig\/$file\n";
 }


# Visulization of MONOD dataset
awk '!/track/ {print $1,$2,$3,$4}' 7-P-1_CpG.bedGraph | sort -k1,1 -k2,2n - > 7-P-1_CpG.bedGraph.Sort.V4
bedGraphToBigWig 7-P-1_CpG.bedGraph.Sort.V4 hg19.chrom.sizes 7-P-1.bw
curl -T 7-P-1.bw ftp://132.239.189.199 --user ucsd002:ucsd002


for i in `ls *bedGraph`
do
awk '!/track/ {print $1,$2,$3,$4}' $i | sort -k1,1 -k2,2n - > $i.S4
done

for i in `ls *.SortV4`
do
bedGraphToBigWig $i hg19.chrom.sizes $i.bw
done

for i in `ls *bw`
do
curl -T $i ftp://132.239.189.199 --user ucsd002:ucsd002
done


foreach my $file (glob("*bw")){
my ($name,undef)=split /\_/,$file;
print "track type=bigWig name=\"$name\" description=\"file\" bigDataUrl=ftp://ucsd002:ucsd002@132.239.189.199/wig/"
}

# utilities tools
bedtools sort -header -i 7-P-1_CpG.bedGraph > 7-P-1_CpG.bedGraph.sort
curl -T 7-P-1.bw ftp://132.239.189.199 --user ucsd002:ucsd002
find mydir -type f -exec curl -u xxx:psw --ftp-create-dirs -T {} ftp://192.168.1.158/public/demon_test/{} \;
track type=bigWig name="My" description="Lab" bigDataUrl=ftp://ucsd002:ucsd002@132.239.189.199/wig/7-P-1.bw
scp 7-P-1.bw shicheng@meangenemachine.dynamic.ucsd.edu:
@132.239.189.199/wig/7-P-1.bw

ftp 132.239.189.199
132.239.189.199/home/ucsd002/wig
samtools tview /home/shg047/oasis/monod/bam/RRBS1/bam/6-P-5.sorted.clipped.bam /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa -p chr17:75283831-75283978
chr17:75283831-75283978 6-P-5.hapInfo.txt
chr17:75283978

setwd("/oasis/tscc/scratch/shg047/monod/mhl")
load("MONOD-Apr6.MHL.RData")
x1<-grep("chr19:58220439-58220459",rownames(data))
x2<-grep("chr19:58220481-58220515",rownames(data))
x3<-grep("chr19:58220626-58220668",rownames(data))
data[x1,]
data[x2,]
data[x3,]

[File:Pancreatic.cancer.plasma.vs.normal.plasma.significant.txt]
[[File:Lung.cancer.plasma.vs.normal.plasma.significant.txt]]
[[File:Colon.cancer.plasma.vs.normal.plasma.significant.txt]]

[[File:13E8.tm-colon.png]]
[[File:146E.tm-lung.png]]
[[File:1508.tm-pancreatic.png]]

SRX381713_normal_lung

less -S monod.mhl.List1.march30.txt

chr10:100027918-100027944

# 1
grep chr10:101302427-101302727 ../hapinfo/6-P-1.hapInfo.txt
grep chr10:100227698-100227747 ../hapinfo/SRX381713_normal_lung.
chr10:100095550-100096429
1E12P20_miPS 0.555555

/home/shg047/oasis/db/mm9/chr10.fa
samtools tview /home/shg047/oasis/db/mm9/chr10.fa

samtools bedcov  /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  /home/shg047/oasis/monod/bam/RRBS2/RRBS2/RRBS-6P16.sorted.clipped.bam 

echo "samtools bedcov /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles/RRBS-6P17.sorted.clipped.bam" | qsub -q glean -N test
How to install Perl package
perl -MCPAN -e "install Getopt::Long"
perl -MCPAN -e "use Getopt::Long"

perl GeneSearch.pl -i pancreatic.ncbi.txt -g GeneSymbolList.txt
/home/shg047/perl5/perlbrew/Getopt-Long-2.48/lib/Getopt/Long.pm
 
bedtools intersect -wao -a colon.cancer.plasma.vs.normal.plasma.significant.txt -b hg19_refGene.bed > colon.cancer.sig.mhl.sorted.Annotated.txt
awk '{print $1}' most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.txt > most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed
perl -p -i -e "s/[:-]/\t/g" most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed 
bedtools intersect -wao -a most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed -b hg19_refGene.bed > most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed.txt.Annotated.txt
awk '{print $8}' most.signficant.66.colon.cancer.plasma.vs.normal.plasma.significant.bed.txt.Annotated.txt > most.signficant.Colon.Gene.list.txt

[[758D.tm.png|400px]]
[[1170.tm-lung-cancer-plasma.png|400px]]
[[3D4E.tm.-pancreatic-cancer-plasma.png|400px]]

chr10:123-456
samtools tview 
/home/shg047/oasis/DennisLo2015/sortbam/BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.PileOMeth.MF.txt

file="lung.cancer.plasma.vs.normal.plasma.significant.txt"
awk '{print $1}' $file > $file.bed
perl -p -i -e "s/[:-]/\t/g" $file.bed
bedtools intersect -wao -a $file.bed -b hg19_refGene.bed > $file.Annotated.txt
awk '{print $8}' $file.Annotated.txt | sort -u > $file.Gene.list.txt

d1<-read.table("xx.lung.target.txt",sep="\t")
d2<-read.table("lung.cancer.plasma.vs.normal.plasma.significant.txt.Gene.list.txt",sep="\t")
d2[na.omit(match(d1[,1],d2[,1])),]

file="pancreatic.cancer.plasma.vs.normal.plasma.significant.txt"
awk 'NR !=1 {print $1}' $file > $file.bed
perl -p -i -e "s/[:-]/\t/g" $file.bed
bedtools intersect -wao -a $file.bed -b hg19_refGene.bed > $file.Annotated.txt
awk '{print $8}' $file.Annotated.txt | sort -u > $file.Gene.list.txt

cd /home/shg047/oasis/DennisLo2015/mhl
awk '{print $1}' dennis.hapinfo2mhl.march28.txt > Dennis.bed
perl -p -i -e "s/[:-]/\t/g" Dennis.bed
/home/shg047/oasis/DennisLo2015/mhl/Dennis.bed

bedtools intersect -wb 

#!/usr/bin/perl
use Cwd;
use strict;
my $file="/home/shg047/oasis/DennisLo2015/mhl/Dennis.bed";
open F,$file;
while(<F>){
chomp;
system("PileOMeth extract -p 5 -q 10 --minDepth 1 -r $_ /home/shg047/oasis/db/hg19/hg19.fa BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam -o Tmp");
my $amf=qx/grep -v track Tmp_CpG.bedGraph | awk '{cov+=$5+\$6;nC+=\$5}END{print nC\/cov}'/;
print "$_\t$amf";
}

less -S ../mhl/dennis.hapinfo2mhl.march28.txt
PileOMeth extract -p 5 -q 10 --minDepth 1 -r chr10:120967293-120967420 /home/shg047/oasis/db/hg19/hg19.fa BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam -o Tmp2
grep -v track Tmp_CpG.bedGraph | awk '{cov+=$5+$6;nC+=$5}END{print nC/cov}' 

PileOMeth extract -p 5 -q 10 --minDepth 1 -r chr10:100069336-100069468 /home/shg047/oasis/db/hg19/hg19.fa BMT1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam -o Tmp
grep -v track Tmp_CpG.bedGraph | awk '{cov+=$5+$6;nC+=$5}END{print nC/cov}' 


chr10:3480-6970
chr10:12345-12347

PileOMeth extract -r chr10:123-456 genome.fa alignments.bam
PileOMeth extract -l a.txt genome.fa alignments.bam  >  c.txt


cd /home/shg047/oasis/monod/bam/RRBS1
perl PileOMethPBS.pl
cd /home/shg047/oasis/monod/bam/RRBS2
perl PileOMethPBS.pl

PileOMeth extract --minDepth 10 /home/shg047/oasis/db/hg19/hg19.fa /home/shg047/oasis/monod/bam/RRBS1/bam/PC-T-7.sorted.clipped.bam -o /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/PC-T-7.sorted.clipped
PileOMeth extract --counts -p 5 -q 10 /home/shg047/oasis/db/hg19/hg19.fa /home/shg047/oasis/monod/bam/RRBS1/bam/PC-T-7.sorted.clipped.bam -o /home/shg047/oasis/monod/bam/RRBS1/MF_PileOMeth/PC-T-7.sorted.clipped


wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
gunzip chr.fa.gz
samtools view -bh RRBS-6P27.sorted.clipped.bam chr1:723205-726280 > Hanne.test.bam
samtools sort -o Hanne.test.sort.bam Hanne.test.bam
samtools index Hanne.test.sort.bam
PileOMeth extract -r chr1:723205-725116  chr1.fa  RRBS-6P27.sorted.clipped.bam  -o y1.ouput.txt
PileOMeth extract -r chr1:725240-726280  chr1.fa  Hanne.test.sort.bam  -o y2.ouput.txt
PileOMeth extract -l interval.input.txt chr1.fa  Hanne.test.sort.bam  -o y1y2.ouput.txt



grep -v track PileOMeth.output | awk '{cov+=$5+$6;nC+=$5}END{print nC/cov}'

samtools tview /home/shg047/oasis/monod/bam/RRBS1/bam/NC-P-12.sorted.clipped.bam /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa -p chr1:725240-725280

NC-P-12.sorqted.clipped_CpG.bedGraph

chr1    725240  725241  0       0       35
chr1    725241  725242  100     1       0
chr1    725255  725256  0       0       36
chr1    725256  725257  100     1       0
chr1    725269  725270  0       0       45
chr1    725270  725271  100     1       0
chr1    725274  725275  4       2       47
chr1    725275  725276  100     1       0
chr1    725279  725280  5       1       16

Rscript --vanilla ~/bin/matrixPlot.R -i "monond-colon-plasma.txt" -o "monond-colon-plasma"
bed="/home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed"
fa="RRBS-6P11 /home/shg047/oasis/db/hg19/hg19.fa"
PileOMeth extract --fraction -q 10 -p 5 -l $bed -o  RRBS-6P11 $fa RRBS-6P11.sorted.clipped.bam  
PileOMeth extract -l /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  RRBS-6P11 /home/shg047/oasis/db/hg19/hg19.fa RRBS-7P30.sorted.clipped.bam
PileOMeth extract -l /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed /home/shg047/oasis/db/hg19/hg19.fa RRBS-7P30.sorted.clipped.bam

 
head /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed > xx.bed



PileOMeth extract /home/shg047/oasis/db/hg19/hg19.fa RRBS-7P30.sorted.clipped.bam


cp /home/shg047/oasis/Tumor-WGBS/HCT116-SRX669642/mergeHapinfo/* ~/oasis/monod/hapinfo
cp /home/shg047/oasis/Tumor-WGBS/HCC-SRX332736/mergeHapinfo/* ~/oasis/monod/hapinfo
cd ~/oasis/monod/hapinfo


perl ~/bin/samInfoPrep4Bam2Hapinfo.pl  .sorted.clipped.bam > ../Saminfo4bam2hapinfo.txt
perl ~/bin/bam2hapInfo2PBS.pl ../Saminfo4bam2hapinfo.txt submit nonbismark



 #!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=4:00:00
 #PBS -o hapinfo2mhl.o
 #PBS -e hapinfo2mhl.e
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/SALK/mergeHapinfo
 perl ~/bin/hapinfo2mhl.pl Salk.HapinfoList > Salk.hapinfo2mhl.march29.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Salk.hapinfo2mhl.march28.txt" -o "Salk.hapinfo2mhl.march28"
 perl ~/bin/hapinfo2mf.pl ./ > Salk.hapinfo2mf.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Salk.hapinfo2mf.march28.txt" -o "Salk.hapinfo2mf.march28"
 
  Rscript --vanilla ~/bin/matrixPlot.R -i "monod.mhl.List1.march29.txt" -o "monod.mhl.List1.march29"

 

 
 
head -n 50 /home/shg047/oasis/AgeBlood/hapinfo/WB-centenarian.hapInfo.txt
head -n 50 /home/shg047/oasis/monod/hapinfo/WB-centenarian.hapInfo.txt

#!/usr/bin/perl
use strict;
use Cwd;

my $dir1="/home/shg047/oasis/Estellar2016/mergeHapinfo/";
chdir $dir1;
my @hap1=glob("$dir1/*lung*hapInfo.txt");
my @hap2=glob("$dir1/*colon*hapInfo.txt");

my $dir2="/home/shg047/oasis/DennisLo2015/hapinfo"
chdir $dir2;
my @hap3=glob("CTR*hapInfo.txt");
my @hap4=glob("Pregnancy*hapInfo.txt");

my $dir3="/home/shg047/oasis/monod/hapinfo";
chdir $dir3;
my @hap5=glob("*hapInfo.txt");

my @file=(@hap1,@hap2,@hap3,@hap4,@hap5);
foreach my $file(@file){
print "$file\n";
}

 
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mhl.march28.txt" -o "Estellar2016.hapinfo2mhl.march28"

 



LC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381719_squamous_cell_tumor_lung.hapInfo.txt
LC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381722_small_cell_tumor_lung.hapInfo.txt
LC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381716_adenocarcinoma_lung.hapInfo.txt
LN	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381713_normal_lung.hapInfo.txt
NC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381553_normal_colon.hapInfo.txt
CC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381569_tumor_colon.hapInfo.txt
CC	/home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381585_metastasis_colon.hapInfo.txt


/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/
/home/shg047/oasis/Estellar2016/mergeHapinfo/



cp /oasis/tscc/scratch/ddiep/Ziller_BAMfiles/Colon_Tumor_Primary* ./

#!/bin/csh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -o wget.log
#PBS -e wget.err
#PBS -M diep.hue.dinh@gmail.com
#PBS -m abe
cd /home/ddiep/dinh_working/Tumor_WGBS
#colon tumor primary tissue
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR949/SRR949210/SRR949210.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR949/SRR949211/SRR949211.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR949/SRR949212/SRR949212.sra
#hct116
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR153/SRR1536575/SRR1536575.sra

GSE46644	GSM1204465	SRX332736	Colon primary tumor
GSE60106	GSM1465024	SRX669642	HCT116

GSE16256
GSE17312
GSE31971
GSE30340


Ziller et al paper.CD14 CD56 CD19

#!/bin/csh
#PBS -l nodes=1:ppn=1
#PBS -q hotel
#PBS -l walltime=12:00:00
#PBS -o wget.log
#PBS -e wget.err
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
# Colon_Tumor_Primary:HCC-SRX332736
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949210/SRR949210_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949210/SRR949210_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949211/SRR949211_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949211/SRR949211_2.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949212/SRR949212_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949212/SRR949212_2.fastq.gz &

# HCT116-SRX669642
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/005/SRR1536575/SRR1536575_1.fastq.gz &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/005/SRR1536575/SRR1536575_2.fastq.gz &

cp /home/shg047/oasis/SALK/mergeHapinfo/STL*.txt ./
wc -l STL003AO-01.hapInfo.txt

# 2016-03-28
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279516/GSM1279516_CpGcontext.Brain_W.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279517/GSM1279517_CpGcontext.Breast.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279518/GSM1279518_CpGcontext.CD19.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279519/GSM1279519_CpGcontext.Colon.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279520/GSM1279520_CpGcontext.Colon_M.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279521/GSM1279521_CpGcontext.Colon_P.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279522/GSM1279522_CpGcontext.H1437.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279523/GSM1279523_CpGcontext.H157.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279524/GSM1279524_CpGcontext.H1672.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279527/GSM1279527_CpGcontext.Lung.txt.gz
ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM1279nnn/GSM1279532/GSM1279532_CpGcontext.U87MG.txt.gz


gzip -frtv9 *

/home/shg047/oasis/Estellar2016/mergeBam


GSE52271


library("GEOquery")
GEOSet <- getGEO("GSE56851")
data <- as.data.frame(exprs(GEOSet[[1]]))
phen <- pData(phenoData(GEOSet[[1]]))

 #!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=4:00:00
 #PBS -o hapinfo2mhl.o
 #PBS -e hapinfo2mhl.e
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/mhl
 Rscript --vanilla ~/bin/matrixPlot.R -i "monod.mhl.march29-List1.txt" -o "monod.mhl.march29-List1"
 cd /home/shg047/oasis/Estellar2016/hapinfo
 perl ~/bin/hapinfo2mhl.pl ./ > Estellar2016.hapinfo2mhl.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mhl.march28.txt" -o "Estellar2016.hapinfo2mhl.march28"
 perl ~/bin/hapinfo2mf.pl ./ > Estellar2016.hapinfo2mf.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mf.march28.txt" -o "Estellar2016.hapinfo2mf.march28"

 
  
 cd /home/shg047/oasis/monod/mhl
 perl ~/bin/hapinfo2mhl.pl hapinfo2mhl.sampleList > monod.mhl.march29.txt
 
 vim hapinfo2mhl-march29.job
 
 
 cd /home/shg047/oasis/Estellar2016/hapinfo
 perl ~/bin/hapinfo2mhl.pl ./ > Estellar2016.hapinfo2mhl.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mhl.march28.txt" -o "Estellar2016.hapinfo2mhl.march28"
 perl ~/bin/hapinfo2mf.pl ./ > Estellar2016.hapinfo2mf.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.hapinfo2mf.march28.txt" -o "Estellar2016.hapinfo2mf.march28"

 
 
 
cd /home/shg047/oasis/monod/hapinfo/hapinfo
cd /home/shg047/oasis/monod/hapinfo/WGBS





[[File:66F0.tm-mhl-dennislo.png]]
[[File:6681.tm.png]]

wc -l dennis.hapinfo2mf.march28.txt
wc -l dennis.hapinfo2mhl.march28.txt

 d1<-read.table("dennis.hapinfo2mf.march28.txt",head=T,row.names=1)
 d2<-read.table("dennis.hapinfo2mhl.march28.txt",head=T,row.names=1)
 library("ggplot2")
 library(RColorBrewer)
 require(KernSmooth)
 pdf("comp-1.smoothscatter.relationship.pdf")
 par(mfrow=c(5,5),mar=c(1,1,1,1))
 for(i in 1:25){
 g = 11
 n=length(d1[,i])
 my.cols <- rev(brewer.pal(g, "RdYlBu"))
 smoothScatter(d1[,i], d2[,i], nrpoints=0.05*n, colramp=colorRampPalette(my.cols), pch=19, cex=.3, col = "green1",xlab="",ylab="")
 }
 dev.off()
 pdf("comp-2.smoothscatter.relationship.pdf")
 par(mfrow=c(2,2),mar=c(2,2,2,2))
 hist(d1[,1],breaks=10,main="MF",col="blue")
 hist(d2[,1],breaks=10,main="MHL",col="blue")
 plot(density(d1[,1],na.rm=T),main="MF",col="blue",lwd=2)
 plot(density(d2[,1],na.rm=T),main="MHL",col="blue",lwd=2)
 dev.off()




# 2016-03-27
vim ~/bin/samInfoPrep4Bam2Hapinfo.pl 
Data<-cbind(GSE53045NormalPBMC,GSE35069NormalPBMC)
cd /home/shg047/oasis/DennisLo2015/sortbam
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ./ > ../Saminfo4bam2hapinfo.txt
perl ~/bin/bam2hapInfo2PBS.pl ../Saminfo4bam2hapinfo.txt submit bismark

vim ~/bin/samInfoPrep4Bam2Hapinfo.pl
cd /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ./ > /home/shg047/oasis/monod/hapinfo/phase2/Saminfo4bam2hapinfo.txt
less -S /home/shg047/oasis/monod/hapinfo/phase2/Saminfo4bam2hapinfo.txt

cd /home/shg047/oasis/monod/hapinfo/phase2/
perl ~/bin/bam2hapInfo2PBS.pl Saminfo4bam2hapinfo.txt submit nonbismark

vim ~/bin/bam2hapInfo2PBS.pl



# 2016-03-26


mkdir test2
samtools view -h /home/shg047/oasis/DennisLo2015/sortbam/T21.1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam chr10:31891611-31891781 > x.sam
samtools view -bh x.sam > x.bam
samtools sort -o x.sort.bam x.bam
samtools index x.sort.bam

samtools tview x.sort.bam /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa -p chr10:31891611-31891781
/home/shg047/bin/mergedBam2hapInfo.pl /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed /oasis/tscc/scratch/shg047/DennisLo2015/sortbam/test2x.sort.bam bismark > ../hapinfo/x.sort.bam.hapInfo.txt

samtools sort -n -o x.nsort.bam x.bam 
bismark_methylation_extractor --bedGraph --zero_based --comprehensive --cutoff 1  --mbias_off --paired-end x.nsort.bam

PileOMeth extract -q 10 -p 5 --mergeContext ~/oasis/db/hg19/hg19.fa x.sort.bam 


/home/shg047/bin/mergedBam2hapInfo.pl /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed /home/shg047/oasis/DennisLo2015/bam/test/test.bam bismark > test.read1.hapInfo.txt
1, Duplicated fragments/reads removing(BS-Seq suitable. RRBS not suitable)
2, Low sequencing quality reads removing (Phred >= 5)
3, Low mapping quality reads removing (MAPQ >= 10)

cd /home/shg047/software/IGV_2.3.68

samtools view -h T21.1.read1_val_1.fq.gz_bismark_bt2_pe.bam | head -n 30001 > x.sam

cd /home/shg047/oasis/DennisLo2015/bam/test
bismark_methylation_extractor --bedGraph --zero_based --comprehensive --cutoff 1  --mbias_off --paired-end test.bam

samtools sort -o test.sort.bam test.bam 
samtools index test.sort.bam
PileOMeth extract --fraction --mergeContext ~/oasis/db/hg19/hg19.fa test.sort.bam 

head test.sort_CpG.bedGraph
head test.bedGraph.gz.bismark.zero.cov

bedtools sort -i test.sort_CpG.bedGraph > test.sort_CpG.sort.bedGraph
bedtools sort -i test.bedGraph.gz.bismark.zero.cov > test.sort.bedGraph.gz.bismark.zero.cov

wc -l test.sort_CpG.sort.bedGraph
wc -l test.sort.bedGraph.gz.bismark.zero.cov

head test.sort_CpG.sort.bedGraph
head test.sort.bedGraph.gz.bismark.zero.cov

bedtools intersect -wo -a test.sort_CpG.sort.bedGraph -b test.sort.bedGraph.gz.bismark.zero.cov > merge.bed

data<-read.table("merge.bed")
pdf("comp.smoothscatter.relationship.pdf")
par(mfrow=c(3,3))
for(i in 1:1){
g = 11
n=length(data[,i])
my.cols <- rev(brewer.pal(g, "RdYlBu"))
smoothScatter(data[,4], data[,10], nrpoints=0.05*n, colramp=colorRampPalette(my.cols), pch=19, cex=.3, col = "green1",xlab="Methylation frequency",ylab="MHL")
}
dev.off()

 chr1:121485049-121485051
 
 
 samtools tview test.sort.bam /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa -p chr1:121485049-121485051



Indx16_S12.sorted.clipped.bam  chr10:102443911-102443400


#!/usr/bin/perl

my $a="0400";
my $a2=0400;
my $b=16;
my $c=$a-$b;
my $c2=$a2-$b;
my $d=0x1a;

if($a){
	print "$a\t$b\t$c\t$a2\t$c2\t$d\n";
}



# 2016-03-24
1, transfer file to new hard-disk:  cp -r /oasis/* /media/NAS2_volume1/shg047
2, fastx_trimmer -Q33 -f 8 -l 43 -i Temp/P608-Tumor.file1.fastq -o Temp/P608-Tumor.file1.trimmed.fastq

qsub -q hotel SRR1035882_1_val_1.fq.gz.job
qsub -q hotel SRR1035895_1_val_1.fq.gz.job

/media/LTS_60T/Dinh

find /tmp -type f -size +50000k -delete



20G     ./bam
1.0K    ./biomarker
9.3G    ./tcga
166G    ./db
33G     ./meth450
1.9G    ./GEO
630G    ./Estellar2016
2.0G    ./mice
398M    ./wbc
911K    ./alice
733M    ./blue
25M     ./haib
2.9G    ./reprogramming
55M     ./ssc
438M    ./bak
422G    ./monod
44G     ./song2
1015M   ./twin
1.3T    ./

rm ../bam/*bam &



# 2016-03-23
cd 

BMT2.read1_val_1.fq.gz_bismark_bt2_pe.bam.bam2mf.err
CTR101_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR103_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR104_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR85_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
HOT240_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err


bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov.sort
bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.sort.bedGraph.gz.bismark.zero.cov 
awk 'NR>8695750 && NR<8695760' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx.bed
4898107
awk '{print $1,$2,$3}' OFS="\t" CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > yy.bed
bedtools sort -i yy.bed > yy.sort.bed
awk '{a[$1];}' xx1.bed | head
bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx.bed
awk '{print FILENAME, NR, FNR, $0}' file1 file2





awk 'NR<4898107' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx1.bed &
awk 'NR>=4898107' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx2.bed &
bedtools sort -i xx1.bed > x1.sort.bed &
bedtools sort -i xx2.bed > x2.sort.bed &
bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > x3.sort.bed 

awk 'NR<8695756' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx1.bed 
bedtools sort -i xx1.bed > x1.sort.bed 
bedtools sort -i CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > x3.sort.bed 



awk 'NR<4898107' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx1.bed &
awk 'NR>4898107' CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov > xx2.bed &

data<-read.table("CTR150_trimmed.fq.gz_bismark_bt2.bedGraph.gz.bismark.zero.cov")


9898107
wc -l CTR150_trimmed.fq.gz_bismark_bt2.sort.bedGraph.gz.bismark.zero.cov 



BMT2.read1_val_1.fq.gz_bismark_bt2_pe.bam.bam2mf.err
CTR101_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR103_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR104_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
CTR85_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
HOT240_trimmed.fq.gz_bismark_bt2.bam.bam2mf.err
Pregnancy.1.read1_val_1.fq.gz_bismark_bt2_pe.bam.bam2mf.err
T21.2.read1_val_1.fq.gz_bismark_bt2_pe.bam.bam2mf.err



# 2016-03-22
samtools view -bh ../../sortbam/Pregnancy.7.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam chr10:2000000-2030000 -o Pregnancy.7.bam
samtools view -bh ../../bam/Pregnancy.7.read1_val_1.fq.gz_bismark_bt2_pe.bam | head -n 3000 > Pregnancy.7.bismark.bam
bismark_methylation_extractor --bedGraph --zero_based --comprehensive Pregnancy.7.csort.bam
bismark_methylation_extractor --bedGraph --zero_based --comprehensive Pregnancy.7.nsort.bam
bismark_methylation_extractor --bedGraph --zero_based --comprehensive Pregnancy.7.bismark.bam


T21.5.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam

HWI-ST1049:8:1101:1199:2048#0/1

#!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test.o
 #PBS -e test.e
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/DennisLo2015/hapinfo
 perl ~/bin/hapinfo2mhl.pl ./ > dennis.hapinfo2mhl.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "dennis.hapinfo2mhl.march28.txt" -o "dennis.hapinfo2mhl.march28"
 perl ~/bin/hapinfo2mf.pl ./ > dennis.hapinfo2mf.march28.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "dennis.hapinfo2mf.march28.txt" -o "dennis.hapinfo2mf.march28"

 
 cd /home/shg047/oasis/monod/hapinfo/hapinfo
 perl ~/bin/hapinfo2mhl.pl ./ > RRBS.Plasma.Phase2.MHL.txt
 Rscript --vanilla ~/bin/matrixPlot.R -i "RRBS.Plasma.Phase2.MHL.txt" -o "RRBS-Plamsa-Phase2"

 
 
 cd /home/shg047/oasis/DennisLo2015/bam/test
 perl ~/bin/bam2hapInfo2PBS.pl test.sam submit bismark > T21.1.read1_val_1.fq.gz_bismark_bt2_pe.mhl

 
 perl ~/bin/hapinfo2mhl.pl ./ submit bismark > T21.1.read1_val_1.fq.gz_bismark_bt2_pe.mhl
 
 
 
 cd /home/shg047/oasis/DennisLo2015/hapinfo
 perl ~/bin/hapinfo2mf.pl ./ > ../dennis.hapinfo2mf.march24.txt
 
 perl ~/bin/hapinfo2mhl.pl ./ > ../dennis.mhl.march24.txt

 cd /home/shg047/oasis/DennisLo2015/bam/test
bismark_methylation_extractor --bedGraph --zero_based --comprehensive CTR98.bam


# 2016-03-22
cd /home/shg047/oasis/DennisLo2015/bam
qsub -q hotel Pregnancy.4.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel Pregnancy.6.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel Pregnancy.7.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel T21.1.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel T21.2.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel T21.3.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job
qsub -q hotel T21.4.read1_val_1.fq.gz_bismark_bt2_pe.bamsort.job

/home/shg047/oasis/AgeBlood/mhl/WB.mhl.march22.txt

ls -lart *bam2sortbam*err
ls -lart *bam2sortbam*err | wc -l 

cd /home/shg047/oasis/Holger2016/fastq_trim/
qsub -q hotel SRR1035734_1_val_1.fq.gz.job 
qsub -q hotel SRR1035723_1_val_1.fq.gz.job
qsub -q hotel SRR1035725_1_val_1.fq.gz.job
qsub -q hotel SRR1035741_1_val_1.fq.gz.job
qsub -q hotel SRR1035738_1_val_1.fq.gz.job
qsub -q hotel SRR1035740_1_val_1.fq.gz.job

cd /home/shg047/oasis/Holger2016/fastq_trim/
qsub -q hotel SRR1035856_1_val_1.fq.gz.job
qsub -q hotel SRR1035829_1_val_1.fq.gz.job
qsub -q hotel SRR1035727_1_val_1.fq.gz.job
qsub -q hotel SRR1035855_1_val_1.fq.gz.job
qsub -q hotel SRR1035799_1_val_1.fq.gz.job

cd /home/shg047/oasis/Holger2016/fastq_trim/
qsub -q hotel SRR1035726_1_val_1.fq.gz.job
qsub -q hotel SRR1035739_1_val_1.fq.gz.job
qsub -q hotel SRR1035735_1_val_1.fq.gz.job
qsub -q hotel SRR1035736_1_val_1.fq.gz.job
qsub -q hotel SRR1035798_1_val_1.fq.gz.job
qsub -q hotel SRR1035804_1_val_1.fq.gz.job
qsub -q hotel SRR1035827_1_val_1.fq.gz.job
qsub -q hotel SRR1035800_1_val_1.fq.gz.job

# 2016-03-21
echo "bsrate -c ~/oasis/db/hg19/hg19.fa -o xx.bsrate xx.rdup" | qsub -q glean
/oasis/tscc/scratch/ddiep/BAMfiles
Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.mhl.march19.txt" -o "Esterllar2016"
PileOMeth extract -l  ~/oasis/db/hg19/hg19.fa





/home/shg047/oasis/DennisLo2015/mr

LC_ALL=C sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -o xx.sort xx.mr
duplicate-remover -S xx_stat.txt -o xx.rdup xx.sort
bsrate -c hg19 -o xx.bsrate xx.rup
grep ˆchrM xx.mr > xx.mr.chrM
bsrate -N -c chrM.fa -o xx.chrM.bsrate xx.mr.chrM
methcounts -n -c hg19 -o xx.meth xx.mr

LC_ALL=C sort -k 1,1 -k 3,3n -k 2,2n -k 6,6 -o xx.sort.end xx.mr


for clinical research. clear data would be very important.
for research the machnism should be great
 
# 2016-03-18
rmapbs -c hg19 -o ../methBam/ SRR949193_1.fastq.gz SRR949193_2.fastq.gz
perl -lane 'print "chr@F[7]\t@F[8]\t@F[8]\t@F[3]\t@F[4]\t@F[6]\t@F[10]\t@F[-3]\t@F[-2]\t@F[-1]" if ! /chrNA/' TCGA.Meth450.esca.ChAMP.DMS.txt  > DMS.bed
grep -v chrNA DMS.bed > DMS2.bed
bedtools intersect -wo -a weimarch20.txt -b DMS2.bed > batch3.region.txt
 
bzip2 -d filename.bz2
This is cutadapt 1.9.1 with Python 2.6.6
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC SRR949197_1.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
gzip: SRR949197_1.fastq.gz: unexpected end of file
cutadapt: error: In read named 'SRR949197.215195112 D1JR8ACXX130107:2:2203:10919:71639 length=99': length of quality sequence (17) and length of read (99) do not match
Cutadapt terminated with exit signal: '256'.
Terminating Trim Galore run, please check error message(s) to get an idea what went wrong...

wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX332%2FSRX332731/SRR949196/SRR949196.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX332%2FSRX332731/SRR949197/SRR949197.sra
wget -r ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/ 

#!/usr/bin/env Rscript
library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

df = read.table(opt$file, header=TRUE)
num_vars = which(sapply(df, class)=="numeric")
df_out = df[ ,num_vars]
write.table(df_out, file=opt$out, row.names=FALSE)

Rscript --vanilla yasrs.R

wget https://cran.r-project.org/src/contrib/getopt_1.20.0.tar.gz
wget https://cran.rstudio.com/src/contrib/optparse_1.3.2.tar.gz
install.packages("getopt_1.20.0.tar.gz")
install.packages("optparse_1.3.2.tar.gz")


 #!/bin/csh
 #PBS -N test
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=14:00:00
 #PBS -o test
 #PBS -e test
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 # cd /home/shg047/oasis/monod/hapinfo
 # perl ~/bin/hapinfo2mhl.pl ./ > ../../mhl.txt
 # wget -r ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/
 cd /home/shg047/oasis/Estellar2016/mergeHapinfo
 perl ~/bin/hapinfo2mf.pl ./ > Estellar2016.mf.march21.txt
 
 cd /home/shg047/oasis/DennisLo2015/mr
 bsrate -c ~/oasis/db/hg19/hg19.fa -o xx.bsrate xx.rup

 
 
 Rscript --vanilla ~/bin/matrixPlot.R -i "Estellar2016.mhl.march19.txt" -o "Esterllar2016"

 cd /home/shg047/oasis/Ziller2013/fastq
 wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/SRR949196_1.fastq.bz2 &
 wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/SRR949196_2.fastq.bz2 &
 wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/SRR949197_1.fastq.bz2 &
 wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA096/SRA096879/SRX332731/SRR949197_2.fastq.bz2 &

 

 perl ~/bin/hapinfo2mhl.pl ./ > ../MHL.OUTPUT.txt


/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles

RRBS-phase2： /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles


ln -s /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles/ bam




cpg_list="/media/Ext12T/DD_Ext12T/BisRef/bisHg19_plusLambda_BWA/hg19_lambda.cpg.positions.txt"
samtools mpileup -BA -f /home/shg047/oasis/db/hg19/hg19.fa $bam_file > $pileup_file
/home/dinh/scripts/BisReadMapper/src/extractMethyl.pl $cpg_list 33 < $pileup_file > $methylFreq_file



# 2016-03-17
ln -s /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles  bam
ln -s /home/k4zhang/my_oasis_tscc/MONOD/Ecker_Tissue_WGBS/BAMfiles bam


/oasis/tscc/scratch/ddiep/Working


SRR949193.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949194.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949195.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949196.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949197.fastq.download.job
-rw------- 1 shg047 k4zhang-group 537 Mar 17 10:56 SRR949198.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949199.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949200.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949201.fastq.download.job
-rw------- 1 shg047 k4zhang-group 455 Mar 17 10:56 SRR949202.fastq.download.job

SRR949193 SRR949202

for i in `seq 194 202`
do 
qsub SRR949$i\fastq.download.job
done






When performing an alignment one must discriminate between different types of bisulfite-treated DNA libraries. In the first, termed directional libraries, adapters are attached to the DNA fragments such that only the original
top or bottom strands will be sequenced. Alternatively, all four DNA strands that arise through bisulfite treatment and subsequent
PCR amplification can be sequenced with the same frequency in nondirectional libraries.


# 2016-03-16
-rw-r--r-- 1 shg047 k4zhang-group          62 Jan 19 12:22 Lymphoma.run3.bam
-rw-r--r-- 1 shg047 k4zhang-group          62 Jan 15 21:56 Lymphoma.run4.bam
-rw-r--r-- 1 shg047 k4zhang-group          62 Feb 10 09:51 Pregnancy.11.bam


qdel 4537516.tscc-mgr.local  shg047      hotel    Lymphoma.run3.re    --      1     16    --   72:00:00 Q       --
qdel 4537527.tscc-mgr.local  shg047      hotel    Lymphoma.run4.re    --      1     16    --   72:00:00 Q       --
qdel 4537636.tscc-mgr.local  shg047      hotel    Pregnancy.11.rea    --      1     16    --   72:00:00 Q       --

chr10:100027918-100027944	TTTT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TTTC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TTCT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCTT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CTTT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TTCC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCTC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CTTC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCCT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CTCT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CCTT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCCC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CTCC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	TCCC	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CCCT	1	100,027,918,100,027,000,000,000,000,000,000,000
chr10:100027918-100027944	CCCC	1	100,027,918,100,027,000,000,000,000,000,000,000

~



qsub Pregnancy.11.read1.fq.gz.job
qsub Lymphoma.run3.read1.fq.gz.job
qsub Lymphoma.run4.read1.fq.gz.job

-rw-r--r-- 1 shg047 k4zhang-group	   62 Jan 19 12:22 Lymphoma.run3.read1_val_1.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group	   62 Jan 15 21:56 Lymphoma.run4.read1_val_1.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group	   62 Feb 10 09:51 Pregnancy.11.read1_val_1.fq.gz_bismark_bt2_pe.bam



bedtools intersect -wo -a znf154.bed -b WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
chr19   58220404	 58220671	 chr19   58220439	 58220459	 chr19:58218498-58222001,B040,4  20
chr19   58220404	 58220671	 chr19   58220481	 58220515	 chr19:58218498-58222001,B042,5  34
chr19   58220404	 58220671	 chr19   58220626	 58220668	 chr19:58218498-58222001,B047,4  42


#!/usr/bin/perl
use strict;
use Cwd;
my $bam_dir=shift @ARGV;
my $bed_dir="/home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed";

chdir $bam_dir;
my @file=glob("*.bam");
foreach my $file(@file){
my ($sample,undef)=split /\./,$file;
print "$sample\t$bam_dir$file\t$bed_dir\n";
}


# 2016-03-14
/home/shg047/monod/predict/phase2
# submit 5 each time. 
qsub SRR949198fastq.download.job
qsub SRR949199fastq.download.job
qsub SRR949200fastq.download.job
qsub SRR949201fastq.download.job
qsub SRR949202fastq.download.job

qsub SRR949203fastq.download.job
qsub SRR949204fastq.download.job
qsub SRR949205fastq.download.job
qsub SRR949206fastq.download.job
qsub SRR949207fastq.download.job

qsub SRR949208fastq.download.job
qsub SRR949209fastq.download.job
qsub SRR949210fastq.download.job
qsub SRR949211fastq.download.job
qsub SRR949212fastq.download.job

qsub SRR949213fastq.download.job
qsub SRR949214fastq.download.job
qsub SRR949215fastq.download.job


# creat bisulfite treated human reference: hg19
cp hg19.fa hg19.c2t.fa
perl -p -i -e 's/CG/M/ig' hg19.c2t.fa
perl -p -i -e 's/C/T/ig' hg19.c2t.fa
perl -p -i -e 's/M/CG/ig' hg19.c2t.fa

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/ 
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/annotation/hg19.c2t.fa -p chr10:100027865
 

chr10:101089382-101089519



Sept 2014: these batch of libraries were generated with the dRRBS protocol (MSP I & Taq I digestion)
Bam folder: /media/Ext12T/DD_Ext12T/RRBS_MONOD/140917_dRRBS/BAMfiles
Mapple_bin_hapInfo folder:
Mld_block_hapInfo folder: /home/kunzhang/CpgMIP/MONOD/Data/140917_dRRBS/mld_block_stringent_hapInfo


cd /home/kunzhang/CpgMIP/MONOD/Data/140917_dRRBS/mld_block_stringent_hapInfo
head NC-P-23.mld_blocks_r2-0.5.hapInfo.txt

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/140917_dRRBS/BAMfiles

# index: http://genome-tech.ucsd.edu/LabNotes/index.php/Dinh/Dinh_2014/NOTES/2014-9-22 
samtools tview Indx16_S12.sorted.clipped.bam /home/shg047/annotation/hg19.c2t.fa -p chr10:102440557-102440826
samtools view Indx16_S12.sorted.clipped.bam  chr10:102443911-102443400


# /home/shg047/oasis/DennisLo2015/bam
cd /home/shg047/oasis/DennisLo2015/fastq/
qsub Lymphoma.run2.read1.fq.gz.job
qsub CTR147.fq.gz.job

cd /home/shg047/oasis/Ziller2013/fastq
qsub SRR949197*job
qsub SRR949196*job
qsub
qsub




# 2016-03-14
File:199.tm.png
File:DF.tm.png

scp /home/shg047/oasis/haib/mhl
cp /home/shg047/monod/hapinfo/STL*.hapInfo.txt ./

setwd("")
/home/shg047/monod/mixHap/hapinfo/mMHL.whole.txt


# 2016-03-11
http://genome-tech.ucsd.edu/LabNotes/index.php/Kun:LabNotes/MONOD/2015-7-6
cd /home/shg047/oasis/monod/mixHap
cut -f 1 colon.data.plsma-hypoall.txt > colon.hypermhl.plasma.txt
cut -f 1 lung.data.plsma.hypoall.txt > lung.hypermhl.plasma.txt
cut -f 1 pancreatic.data.plsma.hypoall.txt > pancreatic.hypermhl.plasma.txt

File:Colon.hypermhl.plasma.txt
File:Lung.hypermhl.plasma.txt
File:Pancreatic.hypermhl.plasma.txt

# 2016-03-10
cd /home/kunzhang/CpgMIP/MONOD/Data/1407-combined_RRBS/mld_blocks_stringent_hapInfo
grep 8054556 7-P-2.mld_blocks_r2-0.5.hapInfo.txt
chr17:8054556-8054584   CCCCC   1	8054556,8054568,8054579,8054581,8054584

cd /home/shg047/monod/hapinfo
grep 8054556 7-P-14.hapInfo.txt
chr17:8054556-8054584   CCCCC   1	8054556,8054568,8054579,8054581,8054584

grep 8054556 7-P-22.hapInfo.txt
chr17:8054556-8054579   CCC     1	8054556,8054568,8054579

grep 95947328 7-P-11.hapInfo.txt
chr9:95947328-95947356  CCCCC   2	95947328,95947337,95947340,95947345,95947356
chr9:95947328-95947356  TTTTT   1	95947328,95947337,95947340,95947345,95947356

grep 100204208 PC-T-2.hapInfo.txt
chr14:100204208-100204257	CCT     1	100204208,100204221,100204232

grep 100204208 PC-P-2.hapInfo.txt
chr14:100204208-100204257	CT      1	100204250,100204257





cd /home/shg047/monod/methyblock/HM450K
bedtools intersect -wa -u -a mhb450bed.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l   # 1551
bedtools intersect -wa -u -a mhb450bed.bed -b mhbbed.bed > mhb450GWBS.bed	      # 1258
bedtools intersect -wa -u -a mhb450GWBS.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l  # 1045

cd /home/shg047/monod/methyblock/encode_rrbs_bed
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l   # 1551
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b /home/shg047/monod/methyblock/HM450K/mhbbed.bed | wc -l    # 8920  
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b /home/shg047/monod/methyblock/HM450K/mhbbed.bed > RRBS_GWBS.bed  
bedtools intersect -wa -u -a RRBS_GWBS.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l  # 7968

cd /home/shg047/monod/methyblock/HM450K
bedtools intersect -wa -u -a mhbbed.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l   # 79704
bedtools intersect -wa -u -a mhbbed.bed -b ~/oasis/db/hg19/CpG.Shore.hg19.bed  | wc -l  # 26103
bedtools intersect -wa -u -a mhbbed.bed -b ~/oasis/db/hg19/CpG. | wc -l   # 3246
bedtools intersect -wa -v -a mhbbed.bed -b ~/oasis/db/hg19/CpGI.hg19.bed | wc -l   # 3246

# 2016-03-09
/home/shg047/oasis/monod/hapinfo/WGBS



CRC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381569_tumor_colon.hapInfo.txt /home/shg047/monod/heatmap
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381716_adenocarcinoma_lung.hapInfo.txt /home/shg047/monod/heatmap
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381719_squamous_cell_tumor_lung.hapInfo.txt /home/shg047/monod/heatmap
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381722_small_cell_tumor_lung.hapInfo.txt /home/shg047/monod/heatmap

2289

CRC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381569_tumor_colon.hapInfo.txt
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381716_adenocarcinoma_lung.hapInfo.txt
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381719_squamous_cell_tumor_lung.hapInfo.txt
LC-T-E001 /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381722_small_cell_tumor_lung.hapInfo.txt

cp /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381569_tumor_colon.hapInfo.txt ./
cp /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381716_adenocarcinoma_lung.hapInfo.txt ./
cp /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381719_squamous_cell_tumor_lung.hapInfo.txt ./
cp /home/shg047/oasis/Estellar2016/mergeHapinfo/SRX381722_small_cell_tumor_lung.hapInfo.txt ./

scp Heatmap.MHL.txt shg047@genome-miner.ucsd.edu:/home/shg047/monod/heatmap/
SRX381569-tumor-colon	CCT
SRX381716-adenocarcinoma-lung	LC
SRX381719-squamous-cell-tumor-lung	LC
SRX381722-small-cell-tumor-lung	LC





# 2016-03-05
mv  *mhl.in.plsma.bed /media/LTS_33T/SG_LTS33T/monod/mixHap/
source("http://bioconductor.org/biocLite.R")
biocLite("DMRcate")
sudo apt-get install iotop
sudo iotop
sudo scp shg047@tscc-login.sdsc.edu:/home/shg047/oasis/monod/hapinfo/SRA/* /home/ucsd002/monod/hapinfo

#!/usr/bin/perl
use strict;
my @file=glob("6-*.hapInfo.txt");
my $i;
foreach my $file(@file){
	$i++;
	my ($cancer,$type,$id)=split /[.-]/g,$file;
	my $id = sprintf("%03d",$i);
	print "$zero_num\n";
	# system("cp $file CRC.$type.$id");
	print "$id\t$i\n";
}


# 2016-03-05
cd /home/shg047/monod/mixHap/mhl
N37-Lung	 Lung
STL001LG-01     Lung
STL002LG-01     Lung
STL002PA-01     Pancreas
N37-Pancreas    Pancreas
STL003PA-01     Pancreas
N37-Colon	Colon
STL001SG-01     Colon
STL003SG-01     Colon


scp shg047@genome-miner.edu:/home/shg047/monod/rrbs_kun/*hapInfo.txt ./
cd /home/shg047/oasis/monod/mhb/hapinfo
cp Colon_primary_tumor.all_chrs.hapInfo.txt HCT116.all_chrs.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/WGBS
cd  /home/shg047/oasis/monod/hapinfo/WGBS 
cd /home/shg047/oasis/monod/haplo/wbc/mergeHapinfo
mv *hapInfo.txt /home/shg047/oasis/monod/haplo/Merge
cd /home/shg047/oasis/monod/haplo/n37/mergeHapinfo
cd /home/shg047/oasis/monod/haplo/salk/mergeHapinfo
cd /home/shg047/oasis/monod/haplo/hesc/mergeHapinfo

11092 Heyn2016.R0.3.methyblock.bed
9296 Heyn2016.R0.4.methyblock.bed
6338 Heyn2016.R0.5.methyblock.bed
3082 Heyn2016.R0.6.methyblock.bed
1366 Heyn2016.R0.7.methyblock.bed
 
 # 2016-03-05
 
 #!/bin/csh
 #PBS -N hapinfo2mhl
 #PBS -q hotel
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o hapinfo2mhl.log
 #PBS -e hapinfo2mhl.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd $PBS_O_WORKDIR
 Rscript MH450DMRSmoking.R
 
 
 perl ~/bin/hapinfo2mhl.pl


# 2016-03-04
Rscript PancancerGSI.R "PancancerMethMatrixaa" "header.txt" "PancancerMethSaminfo_March2016.txt" "PancancerMethMatrixaa.txt"

cd /media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles
cd /oasis/tscc/scratch/ddiep/BAMfiles
perl ~/bin/SaminfoPre4hapinfo.pl > ~/oasis/Estellar2016/SaminfoPre4hapinfo.txt
cd /home/shg047/oasis/Estellar2016/hapinfo
perl ~/bin/bam2hapInfo2PBS.pl ../SaminfoPre4hapinfo.txt

#
for i in `ls *bam`
do 
ls *bam | awk -F_ '{print $1}' | sort -u
done

cd /oasis/tscc/scratch/ddiep/BAMfiles
find ./ -name '*.bam' | {
    read firstbam
    samtools view -h "$firstbam"
    while read bam; do
	 samtools view "$bam"
    done
} | samtools view -ubS - | samtools sort -o merged - 
samtools index merged.bam
ls -l merged.bam merged.bam.bai

-rwx------ 1 shg047 k4zhang-group 1631 Mar  4 09:07 SRX381553_normal_colon.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1601 Mar  4 09:07 SRX381569_tumor_colon.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1751 Mar  4 09:07 SRX381585_metastasis_colon.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1691 Mar  4 09:07 SRX381601_tumor_prostate.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1781 Mar  4 09:07 SRX381611_metastasis_breast.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1631 Mar  4 09:07 SRX381621_tumor_breast.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group 1661 Mar  4 09:07 SRX381631_normal_breast.bam.merge.sh
-rwx------ 1 shg047 k4zhang-group  510 Mar  4 09:07 SRX381646_normal_prostate.bam.merge.sh
ls -larth SRX381646_*
scp shg047@genome-miner.ucsd.edu:/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles/*bam ./
scp shg047@genome-miner.ucsd.edu:/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles/*tumor*bam ./
SRX381646_normal_prostate.bam.merge.sh



cd /media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles

find ./ -name 'SRX381601_tumor_prostate*.bam' | {
read firstbam
echo "$firstbam" 
} 

 samtools view -h ./SRX381601_tumor_prostate.chr14.sorted.clipped.bam

 find .|grep "FooBar"|xargs -I{} cp "{}" ~/foo/bar

 
  ls *bam | grep -v chrLambdaNEB |xargs -I{} cp "{}" /home/shg047/oasis/Estellar2016/bam

  scp shg047@genome-miner.ucsd.edu:/home/shg047/oasis/Estellar2016/bam2/*pl ./
  
 
cd /media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles
find ./ -name 'SRX381725_normal_CD19*.bam' | {
read firstbam
echo "$firstbam"
while read bam; do 
echo "$bam"
done
}

samtools view ./SRX381725_normal_CD19.chr19.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chrX.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr14.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr20.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chrY.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr15.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr7.sorted.clipped.bam | head 
samtools view ./SRX381725_normal_CD19.chrM.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr11.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chrLambdaNEB.sorted.clipped.bam  | head
samtools view ./SRX381725_normal_CD19.chr10.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr9.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr17.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr3.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr21.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr12.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr22.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr18.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr4.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr5.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr13.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr16.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr8.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr2.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr6.sorted.clipped.bam | head
samtools view ./SRX381725_normal_CD19.chr1.sorted.clipped.bam | head




find ./ -name 'SRX381601_tumor_prostate*.bam' | {
read firstbam
echo "$firstbam"
samtools view -h "$firstbam"
while read bam; do
echo "$bam"
samtools view "$bam"
done
 } | samtools view -ubS - | samtools sort - /home/shg047/oasis/Estellar2016/bam/merged
samtools index /home/shg047/oasis/Estellar2016/bam/merged.bam





shg047@tscc-login.sdsc.edu


scp * shg047@tscc-login.sdsc.edu:/home/shg047/oasis/TCGA/Meth/Data


wget --ftp-user='zhang' --ftp-password='HDIIWpP' ftp://169.228.63.66/
wget -r --user='zhang' --password='HDIIWpP' ftp://169.228.63.66/

bedtools intersect -a /home/shg047/db/hg19/encode/encode.*.hg19.bed -b en.mhb.hypo.bed

rm endo.mhb.hypo.tf
rm meso.mhb.hypo.tf
rm ecto.mhb.hypo.tf
for i in `ls /home/shg047/db/hg19/encode/encode.*.hg19.bed`
do
bedtools window -w 100 -a endo.mhb.hypo.bed -b $i >> endo.mhb.hypo.tf
bedtools window -w 100 -a meso.mhb.hypo.bed -b $i >> meso.mhb.hypo.tf
bedtools window -w 100 -a ecto.mhb.hypo.bed -b $i >> ecto.mhb.hypo.tf
done
cat endo.mhb.hypo.tf | awk '{print $7}' | sort -u > endo.mhb.hypo.uni.tf
cat meso.mhb.hypo.tf | awk '{print $7}' | sort -u > meso.mhb.hypo.uni.tf
cat ecto.mhb.hypo.tf | awk '{print $7}' | sort -u > ecto.mhb.hypo.uni.tf


for i in `ls /home/shg047/db/hg19/encode/encode.*.hg19.bed`
do
bedtools intersect -a $i -b ecto.mhb.hypo.bed 
done


/home/shg047/oasis/monod/bin/bam2hapinfo.pl /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.chr6.mld_blocks_r2-0.5.bed  /home/k4zhang/my_oasis_tscc/MONOD/tumor_WGBS/BAMfiles/Colon_primary_tumor.chr7.sorted.clipped.bam > Colon_primary_tumor.chr7.hapInfo.txt

/home/oasis/db/
cd salk
perl ../bam2hapinfoPBS.pl ../Ecker_Tissue_WGBS_sampleInfo_WGBS-pooled-mld-blocks.txt
cd ../tumor_wgbs/
perl ../bam2hapinfoPBS.pl ../tumor_WGBS_sample_info_WGBS-pooled-mld-blocks.txt
cd ../wbc/
perl ../bam2hapinfoPBS.pl ../TSCC_whole_blood_WGBS_sampleInfo_WGBS-pooled-mld-blocks.txt
cd ../n37/
perl ../bam2hapinfoPBS.pl ../N37_10_tissue_sampleInfo_WGBS-pooled-mld-blocks.txt
cd ../hesc/
perl ../bam2hapinfoPBS.pl ../H1ESC_sampleInfo_WGBS-pooled-mld-blocks.txt

HWI-D00506:1:1102:13764:14978#0/3:F     16      chr5    177653119	14      84M     *	0	0	CAAACTAAAATACAATAACGCGATCTCAACTCACTACAACCTCTACCTCTCAAACTCAAACAATTCTCCTACTTCAACCTCCCA    bbbeeeeeggfggiiihhhiifcgfiiihihiiihhi
HWI-D00506:1:1106:21115:64108#0/1_HWI-D00506:1:1106:21115:64108#0/3_121:R	16      chr5    177653119	14      121M    *	0	0	CAAACTAAAATACAATAACTCTATCTCAACTCACTACAACCTCTACCTCTCAAACTCAAACAATTCTCCTACTTCAACCTCCCAA
HWI-D00506:1:1102:13764:14978#0/1:R     16      chr5    177653121	14      84M     *	0	0	AACTAAAATACAATAACGCGATCTCAACTCACTACAACCTCTACCTCTCAAACTCAAACAATTCTCCTACTTCAACCTCCCAAA    bbbeeeeefcggghhgfghiihhiihfhihhigfhii

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/ 
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/annotation/hg19.c2t.fa -p chr5:177653138
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/db/hg19/hg19.fa -p chr5:177653138

# I found MHL at chr5:177653138-177653229 in 6-P-1 is missing, however, there are 3 reads aligned to this region in the bam file. 
# I need check what's wrong with it? Anything should be de-bug for bam2hapinfo.pl? 
mkdir 
samtools view -b /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam chr5:177653138-177653229 -o 6-P-1.test.bam 
perl ../../bin/bam2hapinfo.pl target.bed 6-P-1.test.bam
grep 177653138 6-P-1*
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/db/hg19/hg19.fa -p chr5:177653138
samtools view 6-P-1.test.bam 
samtools view -q 6-P-1.test.bam 


Illumicode  MouseWG-6 v2.0



 #!/bin/csh
 #PBS -N bam2MHB
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=7:00:00
 #PBS -o bam2mhb.log
 #PBS -e bam2mhb.err
 #PBS -V
 #PBS -M shicheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/monod/hapinfo
 perl ../bin/hapinfo2mhl.pl ./




cd /home/shg047/oasis/monod/haplo/n37
perl hapMergeByChrosome.pl
mv N37-CRBL.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Colon.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-FL.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Heart.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Liver.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Lung.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Pancreas.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-SI.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-SM.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv N37-Stomach.hapInfo.txt /home/shg047/oasis/monod/hapinfo/

cd /home/shg047/oasis/monod/haplo/salk
perl ../n37/hapMergeByChrosome.pl ./
mv STL001BL-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001FT-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001GA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001LG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001LV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001PO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001RV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001SB-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001SG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001SX-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL001TH-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002AD-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002AO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002EG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002FT-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002GA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002LG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002OV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002PA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002PO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002SB-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL002SX-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003AD-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003AO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003EG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003FT-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003GA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003LV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003PA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003PO-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003RA-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003RV-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003SB-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003SG-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL003SX-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/
mv STL011LI-01.hapInfo.txt  /home/shg047/oasis/monod/hapinfo/

cd /home/shg047/oasis/monod/haplo/wbc
perl ../n37/hapMergeByChrosome.pl ./
mv centenarian.hapInfo.txt   /home/shg047/oasis/monod/hapinfo/
mv middle-age.hapInfo.txt   /home/shg047/oasis/monod/hapinfo/
mv new-born.hapInfo.txt   /home/shg047/oasis/monod/hapinfo/

cd /home/shg047/oasis/monod/haplo/hesc
perl ../n37/hapMergeByChrosome.pl ./
mv methylC-seq_h1+bmp4_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1+bmp4_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1-msc_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1-msc_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1-npc_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1-npc_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1_mesendoderm_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1_mesendoderm_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1_r1.hapInfo.txt /home/shg047/oasis/monod/hapinfo/
mv methylC-seq_h1_r2.hapInfo.txt /home/shg047/oasis/monod/hapinfo/





cd /home/shg047/monod/methyblock/HM450K/mhb450bed.bed 
bedtools intersect -wa -u -a mhb450bed.bed -b ../WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  # 1258
bedtools shuffle -i /home/shg047/monod/methyblock/HM450K/mhb450bed.bed -incl /home/shg047/oasis/db/hg19/CRGmapability.hg19.exclude.bed -g ~/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/methyblock/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l  # random 48, observe: 
bedtools intersect -wa -u -a /home/shg047/monod/methyblock/HM450K/mhb450bed.bed -b /home/shg047/oasis/monod/methyblock/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l  # random 48, observe: 1258


cd /home/shg047/monod/methyblock/encode_rrbs_bed
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b ../WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l   # 8920
bedtools shuffle -i /home/shg047/monod/methyblock/encode_rrbs_bed/encode.rrbs.mhb.txt -incl /home/shg047/oasis/db/hg19/CRGmapability.hg19.exclude.bed -g ~/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/methyblock/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l  # random 337, observe: 




bedtools shuffle -i /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed -incl /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed -g /home/shg047/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l





/home/shg047/monod/rrbs_kun/*hapInfo.txt

load("Encode.RRBS.MHB.RData")
cor2bed<-function(cor){
  a<-unlist(lapply(strsplit(cor,split=c(":")),function(x) strsplit(x,"-")))
  bed<-matrix(a,ncol=3,byrow=T)
  return(data.frame(bed))
}
rrbsmhb<-cor2bed(rownames(subset1))
write.table(rrbsmhb,file="encode.rrbs.mhb.txt",col.names=F,row.names=F,quote=F,sep="\t")
bedtools intersect -wa -u -a encode.rrbs.mhb.txt -b ../WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
bedtools intersect -wa -u -b encode.rrbs.mhb.txt -a ../WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l

../../mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 6-T-ALL.hapInfo.txt 6-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > RRBS_targets_LOD_CC.txt &
../../mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 7-T-ALL.hapInfo.txt 7-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt  | sort -k 2,2nr  > RRBS_targets_LOD_LC.txt &
../../mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt PC-T-ALL.hapInfo.txt PC-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > RRBS_targets_LOD_PC.txt &

cd /home/kunzhang/CpgMIP/MONOD/Data/1407-combined_RRBS/expanded_analysis_31Aug2014

perl /home/kunzhang/CpgMIP/MONOD/Data/mixMethHapAnalysisStringent_27Aug14.pl NC-P-ALL.hapInfo.txt 6-T-1.hapInfo.txt 6-P-1.hapInfo.txt  /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_RRBS_targets.BED.txt > /home/shg047/monod/mixHap/zhang/6-P-1_tumor_hap.txt

grep 105309703-105309944 NC-P-ALL.hapInfo.txt 
grep 105309703-105309944 6-T-ALL.hapInfo.txt 
grep 105309703-105309944 6-P-ALL.hapInfo.txt




mixMethHapAnalysisStringent_27Aug14

8920/
diff mixMethHapAnalysisStringent_27Aug14.pl mixMethHapAnalysisStringent_02Nov14.pl



my $hap_file_NCP = $ARGV[0];
my $hap_file_Tumor = $ARGV[1];
my $hap_file_PP = $ARGV[2];
my $mC_bed_file_NCP = $ARGV[3];
my $informative_probe_file = $ARGV[4];


my $hap_file_NCP = "NC-P-ALL.hapInfo.txt";
my $hap_file_Tumor = "6-T-ALL.hapInfo.txt";
my $hap_file_PP = "6-P-ALL.hapInfo.txt";
my $mC_bed_file_NCP = "/home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt";
my $informative_probe_file = "";





cd /home/kunzhang/CpgMIP/MONOD/Data/1407-combined_RRBS/expanded_analysis_31Aug2014

perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 6-T-ALL.hapInfo.txt 6-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > /home/kunzhang/CpgMIP/MONOD/DataRRBS_targets_LOD_CC.txt 
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 7-T-ALL.hapInfo.txt 7-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt  | sort -k 2,2nr  > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_LC.txt &
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt PC-T-ALL.hapInfo.txt PC-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_PC.txt &


# 2016-02-26
wget https://cran.r-project.org/src/contrib/mhsmm_0.4.14.tar.gz
install.packages("mhsmm_0.4.14.tar.gz")
source("http://bioconductor.org/biocLite.R")
biocLite("MethylSeekR")
Hidden Markov Model to identify PMDs, UMRs and LMRs

/home/shg047/oasis/monod/bin/LDR2extent.pl chr5:112039140-112047142 APC.chr5.sqr > APC.full.region.sqr
head -n 100 /home/kunzhang/HsGenome/hg19/HsGenome19.CpG.positions.txt > head.txt
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  APC.ext chr5:112043054-112197528 < /home/shg047/oasis/monod/haplo/hapinfo.txt
perl /home/shg047/bin/LDR2extent.pl chr5:112043054-112197528 APC.chr5.rsq /home/kunzhang/HsGenome/hg19/HsGenome19.CpG.positions.txt > APC.chr5.rsq.ext
my $cpg_position_file="/home/shg047/oasis/db/hg19/meth/bismark/HsGenome19.CpG.positions.txt";
# APC regions
cd /home/shg047/oasis/monod/haplo/
cat /home/shg047/oasis/monod/haplo/n37/*hapInfo.txt >> hapinfo.txt
cat /home/shg047/oasis/monod/haplo/wbc/*hapInfo.txt >> hapinfo.txt
cat /home/shg047/oasis/monod/haplo/hesc/*hapInfo.txt >> hapinfo.txt
cat /home/shg047/oasis/monod/haplo/tumor_wgbs/*hapInfo.txt >> hapinfo.txt
cat /home/shg047/oasis/monod/haplo/salk/*hapInfo.txt >> hapinfo.txt

cd /home/shg047/oasis/monod/haplo/
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  APC chr5:111986595-112228720 < /home/shg047/oasis/monod/haplo/hapinfo.txt
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  ZNF154 chr19:58206356-58222715 < /home/shg047/oasis/monod/haplo/hapinfo.txt
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  SHOX2 chr3:157808200-157829413 < /home/shg047/oasis/monod/haplo/hapinfo.txt
perl ~/oasis/monod/bin/haploinfo2LDR2.pl  DIRAS3 chr1:68511645-68516481 < /home/shg047/oasis/monod/haplo/hapinfo.txt


awk 'NR==9467146 {print}' hapinfo.txt
awk 'NR==8693459 {print}' hapinfo.txt
awk 'NR==8693437 {print}' hapinfo.txt

chr5:111986595-112228720
.
# 2016-02-24
SRR1035775
SRR1035774
SRR1035786
SRR1035787
SRR1035788
SRR1035818
SRR1035778
SRR1035893
SRR1035815
SRR1035776
SRR1035893
SRR1035814
SRR1035774
SRR1035777
SRR1035813
SRR1035777
SRR1035814
SRR1035773
SRR1035788
SRR1035787
SRR1035775
SRR1035773
SRR1035777
SRR1035786
SRR1035786
SRR1035776
SRR1035778
SRR1035822
SRR1035815

􀀛􀀚􀀘􀀃􀀫􀁒􀁗􀁈􀁏􀀃􀀦􀁌􀁕􀁆􀁏􀁈􀀃􀀶􀁒􀁘􀁗􀁋􀀏􀀃􀀰􀁌􀁖􀁖􀁌􀁒􀁑􀀃􀀹􀁄􀁏􀁏􀁈􀁜􀀏􀀃􀀶􀁄􀁑􀀃􀀧􀁌􀁈􀁊􀁒


bismark_methylation_extractor -p --no_overlap --ignore_r2 2 --comprehensive --gzip --report --multicore 4 --bedGraph --ample_memory --cytosine_report --CX --genome_folder ~/genomes/Heinz_pennellii_organelles_F1_genome/ --split_by_chromosome -o methyl_extraction/ P1.deduplicated.bam

# cytosineReport to MethylKit
awk  '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1,$2,$3,$4/($4+$5),$4+$5;}' cytosineReport.txt > outputForMethylKit

fastq-dump --split-files --skip- --gzip $id 

for i in {4391177..4391187}; do qdel $i.tscc-mgr.local; done

# SRR390728 is pair-end
fastq-dump -X 20 -Z SRR390728
fastq-dump --split-files  -X 20 -Z SRR390728 
fastq-dump --split-3  -X 20 -Z SRR390728 

fastq-dump -X 20 SRR390728
fastq-dump --split-files  -X 20 SRR390728 
fastq-dump --split-3  -X 20 SRR390728 
fastq-dump -X 5 -Z SRR390728

# SRR2542443 is single-end
fastq-dump --split-files -X 20 SRR2542443 

cd /home/shg047/oasis/Holger2016
screen -S SRA112056.Download 
wget -r -c -l 2 -nH --cut-dirs=4 ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA112/SRA112056/
qwget -r -c -l 2  -nH  ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA112/SRA112056/
wget -r -l 2 ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA112/SRA112056

# 2016-02-21
scp shg047@genome-miner.ucsd.edu:/home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt ./
cd /home/kunzhang/CpgMIP/MONOD/Data/
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 6-T-ALL.hapInfo.txt 6-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_CC.txt 
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 7-T-ALL.hapInfo.txt 7-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt  | sort -k 2,2nr  > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_LC.txt &
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt PC-T-ALL.hapInfo.txt PC-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > /home/shg047/monod/mixHap/zhang/RRBS_targets_LOD_PC.txt &

/home/kunzhang/CpgMIP/MONOD/Data/mixMethHapAnalysis_27Oct14.pl
/home/kunzhang/CpgMIP/MONOD/Data/mixMethHapAnalysisStringent_02Nov14.pl

cd /home/kunzhang/CpgMIP/MONOD/Data/1407-combined_RRBS/expanded_analysis

scp shg047@genome-miner.ucsd.edu:/home/kunzhang/CpgMIP/MONOD/Data/*pl ./
scp shg047@genome-miner.ucsd.edu:/home/kunzhang/CpgMIP/*pl ./

cd /home/kunzhang/CpgMIP


mixMethHapAnalysisStringent_27Aug14.pl

cd /home/shg047/monod/hapinfo
cat /home/shg047/monod/rrbs_kun/6-P-*.hapInfo.txt > ../mixHap/6-P-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/7-P-*.hapInfo.txt > ../mixHap/7-P-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/PC-P-*.hapInfo.txt > ../mixHap/PC-P-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/NC-P-*.hapInfo.txt > ../mixHap/NC-P-ALL.hapInfo.txt

cat /home/shg047/monod/rrbs_kun/6-T-*.hapInfo.txt > ../mixHap/6-T-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/7-T-*.hapInfo.txt > ../mixHap/7-T-ALL.hapInfo.txt
cat /home/shg047/monod/rrbs_kun/PC-T-*.hapInfo.txt > ../mixHap/PC-T-ALL.hapInfo.txt

cd /home/shg047/monod/mixHap
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 6-T-ALL.hapInfo.txt 6-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > RRBS_target_LOD_CC.txt 
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt 7-T-ALL.hapInfo.txt 7-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt  | sort -k 2,2nr  > RRBS_target_LOD_LC.txt 
perl /home/shg047/monod/bin/mixMethHapAnalysis_19Aug14.pl  NC-P-ALL.hapInfo.txt PC-T-ALL.hapInfo.txt PC-P-ALL.hapInfo.txt /home/kunzhang/CpgMIP/MONOD/Public_data/WB_middle_age_UMR.BED.txt | sort -k 2,2nr > RRBS_target_LOD_PC.txt 
./find_NC-P_high-LOD_targets.pl > NC-P_high-LOD_targets.txt
  
  
find /home/kunzhang/CpgMIP/MONOD/Data/ -name *pl | grep find_NC-P_high
cp /home/kunzhang/CpgMIP/MONOD/Data/*.sh ./


# creat bisulfite treated human reference: hg19
cp hg19.fa hg19.c2t.fa
perl -p -i -e 's/CG/M/ig' hg19.c2t.fa
perl -p -i -e 's/C/T/ig' hg19.c2t.fa
perl -p -i -e 's/M/CG/ig' hg19.c2t.fa

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/ 
samtools tview  /media/Ext12T/DD_Ext12T/RRBS_MONOD/Bam_Merged/6-P-1.merged.bam /home/shg047/annotation/hg19.c2t.fa -p chr10:100027865

chr10:101089382-101089519
Sept 2014: these batch of libraries were generated with the dRRBS protocol (MSP I & Taq I digestion)
Bam folder: /media/Ext12T/DD_Ext12T/RRBS_MONOD/140917_dRRBS/BAMfiles
Mapple_bin_hapInfo folder:
Mld_block_hapInfo folder: /home/kunzhang/CpgMIP/MONOD/Data/140917_dRRBS/mld_block_stringent_hapInfo


cd /home/kunzhang/CpgMIP/MONOD/Data/140917_dRRBS/mld_block_stringent_hapInfo
head NC-P-23.mld_blocks_r2-0.5.hapInfo.txt

cd /media/Ext12T/DD_Ext12T/RRBS_MONOD/140917_dRRBS/BAMfiles

# index: http://genome-tech.ucsd.edu/LabNotes/index.php/Dinh/Dinh_2014/NOTES/2014-9-22 
samtools tview Indx16_S12.sorted.clipped.bam /home/shg047/annotation/hg19.c2t.fa -p chr10:102440557-102440826
samtools view Indx16_S12.sorted.clipped.bam  chr10:102443911-102443400


chr10:101281221-101281239
chr10:101281274-101281291
chr10:101281879-101281896 
# 2016-02-19
head /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed 

/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge/BAMfiles

cpg_list="/media/Ext12T/DD_Ext12T/BisRef/bisHg19_plusLambda_BWA/hg19_lambda.cpg.positions.txt"
samtools mpileup -BA -f hg19.fa $bam_file > $pileup_file
/home/dinh/scripts/BisReadMapper/src/extractMethyl.pl $cpg_list 33 < $pileup_file > $methylFreq_file

awk '{print $3-$2}' /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed

awk '{print $3-$2}' Haib.merge_RD10_80up.mld_blocks_r2-0.1.bed | sort -u
awk '{print $3-$2}' /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  |sort -n -u
awk '{print $3-$2}' /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed  | sort -n -u

scp shg047@genome-miner.ucsd.edu:/home/shg047/monod/methyblock/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed ./
scp shg047@genome-miner.ucsd.edu:/home/shg047/monod/methyblock/*.bed ./

bedtools intersect -wa -u -a  /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed  | wc -l

# test 1
wc -l /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed
bedtools intersect -wa -u -a /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l 
wc -l /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed
bedtools shuffle -i /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed -incl /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed -g /home/shg047/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed | wc -l

# test 2
wc -l /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed
bedtools intersect -wa -u -a /home/shg047/oasis/Haib/mhb/haib.RD10_80up.genomecov.bed -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l 
wc -l /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed
bedtools shuffle -i /home/shg047/oasis/Haib/mhb/Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed -excl -g /home/shg047/oasis/db/hg19/hg19.chrom.sizes | bedtools intersect -wa -u -a - -b /home/shg047/oasis/monod/mhb/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed | wc -l


# Cancer specific methylation haplotype
http://genome-tech.ucsd.edu/LabNotes/index.php/Cancer_specific_methylation_haplotype
http://genome-tech.ucsd.edu/LabNotes/index.php/Shicheng:Calendar/NOTES/2016-2-1
http://genome-tech.ucsd.edu/LabNotes/index.php/Kun:LabNotes/MONOD/2014-10-8_WGBS_RRBS
http://genome-tech.ucsd.edu/LabNotes/index.php/Kun:LabNotes/MONOD/2014-10-8_WGBS_RRBS


less /home/shg047/oasis/Holger2016/sra/SRR1035896.trim.err
 scp shg047@genome-miner.ucsd.edu:/home/kunzhang/CpgMIP/MONOD/Data/methHapClassfier.pl ./
 methHapClassfier.pl  N37-Cerebellum.chr1.hapInfo.txt > N37-Cerebellum.methHapCounts.txt
 cat PC-T-*.methHapCounts.txt | ../scripts/get_tumor_specific_HMH_regions.pl NC-P-ALL_RRBS-dRRBS.methHapCounts.txt > PC-T_NC-plasma_specific_HMH_regions.bed
 cat 7-T-*.methHapCounts.txt | ../scripts/get_tumor_specific_HMH_regions.pl NC-P-ALL_RRBS-dRRBS.methHapCounts.txt > LC-T_NC-plasma_specific_HMH_regions.bed
 cat 6-T-*.methHapCounts.txt CTT-frozen-100ng.methHapCounts.txt  | ../scripts/get_tumor_specific_HMH_regions.pl NC-P-ALL_RRBS-dRRBS.methHapCounts.txt > CRC-T_NC-plasma_specific_HMH_regions.bed

Question 1: Cannot locate ... in @INC - Perl Maven
Question 2: How to install the module
Question 3: Where I have installed my module 
Question 4: how to load module

Answer 1:
module path is not in the @INC. You need add the path to @INC

Answer 2:  
cpan
install Sort::Array

Answer 3:  
perldoc -l XML::Simple
perldoc -l Sort::Array

Answer 4:
export PERL5LIB=$PERL5LIB:/home/shg047/perl5/perlbrew/perls/perl-5.22.0/lib/site_perl/5.22.0/Sort/
export PERLLIB=$PERLLIB:/home/shg047/perl5/perlbrew/perls/perl-5.22.0/lib/site_perl/5.22.0/Sort/
source ~/.bashrc 


qsub SRR1035860.job
qsub SRR1035849.job
qsub SRR1035856.job
qsub SRR1035860.job
qsub SRR1035863.job
qsub SRR1035869.job
qsub SRR1035879.job
qsub SRR1035888.job
qsub SRR1035895.job

 #!/bin/csh
 #PBS -N bam2MHB
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=72:00:00
 #PBS -o bam2mhb.log
 #PBS -e bam2mhb.err
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/Haib/sortBam
 ../mergedBam2hapInfo.pl ./haib.RD10_80up.genomecov.bed haib.merge.sort.bam > Haib.merge.RD10_80up.hapinfo.txt  # get hapinfo
 ../hapInfo2mld_block.pl ./Haib.merge.RD10_80up.hapinfo.txt 0.5 >  Haib.merge_RD10_80up.mld_blocks_r2-0.5.bed   # get methylation block


 #!/bin/csh
 #PBS -n CTR97_trimmed.fq.gz_bismark_bt2.sort.bam
 #PBS -q pdafm
 #PBS -l nodes=1:ppn=16
 #PBS -l walltime=72:00:00
 #PBS -o CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.log
 #PBS -e CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.err
 #PBS -V
 #PBS -M shihcheng.guo@gmail.com
 #PBS -m abe
 #PBS -A k4zhang-group
 cd /home/shg047/oasis/Haib/sortBam
 # samtools cat -h header.sam -o haib.merge.bam *sort.bam
 samtools sort -@ 16 haib.encode.merge.bam -o haib.merge.sort.bam
 samtools index haib.merge.sort.bam
 bedtools genomecov -bg -split -ibam haib.merge.sort.bam >   haib.merge.bam.pool.bed
 awk '$4>9 { print $1"\t"$2"\t"$3}'  haib.merge.bam.pool.bed | bedtools merge -d 10 -i - > haib.RD10.genomecov.bed
 awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' haib.RD10.genomecov.bed > haib.RD10_80up.genomecov.bed

 for i in {1..22,"X","Y"} 
 do 
 grep chr1 haib.RD10_80up.genomecov.bed |awk '{sum+=($3-$2)} END { print "chr1 Sum = ", sum, " Average = ",sum/NR, "N = ", NR}' >> N37_WGBS_tumor_seqCap_RRBS_tumor_NC_RD10_80up.mld_blocks.summary.txt
 done
 
  for i in `seq 1 12`;
  do
  grep chr$1 haib.RD10_80up.genomecov.bed |awk '{sum+=($3-$2)} END { print "chr1 Sum = ", sum, " Average = ",sum/NR, "N = ", NR}' >> Haib_RD10_80up.mld_blocks.summary.txt
  done 
 
 
 mkdir /home/shg047/oasis/Haib/mhb
 mv 
 genome.cov.sh					 

mv haib.encode.merge.bam ../mhb			 
mv CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.err  ../mhb   
mv haib.merge.sort.bam	    ../mhb		     
mv haib.merge.sort.bam.bai      ../mhb		      
mv haib.merge.bam.pool.bed	 ../mhb		    
mv haib.RD10.genomecov.bed	    ../mhb		 
mv haib.RD10_80up.genomecov.bed	 ../mhb	      
mv CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.log     ../mhb

 
# 2016-2-16
Dinh methylation collection page @ genome-miner
/media/LTS_60T/Dinh/WGBS_LTS33/Hg19/Estellar_Bellvitge
DNA Methylation by Reduced Representation Bisulfite Seq from ENCODE/HudsonAlpha 

cd /home/kunzhang/CpgMIP/MONOD/Data
wget -r -l 2 ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA112/SRA112056


# 2016-2-13
(head -n 1 $file && tail -n +2 $file | sort -k1,1 -k2,2n | awk '{print $1,$2,$3,$4}' OFS="\t" )  > $sample.sort.bedGraph
bedGraphToBigWig NC-P-2.bedGraph_CpG.sort.bedGraph hg19.chrom.sizes  NC-P-2.bw
fetchChromSizes hg19 > hg19.chrom.sizes
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes

echo {0..9} | xargs -n 2

shicheng@meangenemachine2
ifconfig: 132.239.189.199
sudo passwd: fudan1108
sudo apt-get update
sudo apt-get install vim
sudo apt-get install vsftpd
sudo vim /etc/vsftpd.conf
# Listen=YES
anonymous_enables=NO
local_enable=YES
write_enable=YES
xferlog_file=/var/log/vsftpd.log
ftpd_banner= Welcome to my new FTP Serve at UCSD
save and exit
sudo cat /var/log/vsftpd.log
sudo service vsftpd restart
sudo adduser ucsd002  passwd: ucsd002

 
sudo telnet localhost 21
ps -aux | grep vsftpd
sudo service vsftpd restart
sudo netstat -ntaulp | grep vsftpd
sudo adduser genemean001  passwd: genemean001
sudo adduser ucsd002  passwd: ucsd002


 
# 2016-2-12
cd ~/bin/
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigCorrelate
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigEncode
chmod 700 * 
cd -

wigCorrelate

#!/bin/csh
#PBS -n CTR97_trimmed.fq.gz_bismark_bt2.sort.bam
#PBS -q glean
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -o CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.log
#PBS -e CTR97_trimmed.fq.gz_bismark_bt2.sort.bam.err
#PBS -V
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group

cd /home/shg047/oasis/biomark
PileOMeth extract  -q 1 -p 1 --minDepth 1 /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam  -o 7-P-3.PileOMeth.bedGraph
samtools depth -b ~/oasis/db/hg19/meth/bismark/HsGenome19.CpG.positions.txt /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam > 7-P-3.coverage
samtools tview -p chr1:10923 7-P-3.sorted.clipped.bam ~/oasis/db/hg19/meth/bismark/hg19.fa

5697312 7-P-3.coverage

37076

samtools tview -p chr1:10923 7-P-3.sorted.clipped.bam ~/oasis/db/hg19/meth/bismark/hg19.fa


[shg047@tscc-login1 biomark]$ wc -l 7-P-3.PileOMeth.bedGraph_CpG.bedGraph
37076 7-P-3.PileOMeth.bedGraph_CpG.bedGraph

/home/shg047/oasis/monod/bam
less /home/shg047/oasis/monod/bam/*coverage

chr1    10854   15
chr1    10857   16
chr1    10860   23
chr1    10866   23
chr1    10884   24
chr1    10886   24
chr1    10902   24
chr1    10907   24
chr1    10913   24
chr1    10923   23
chr1    10928   23
chr1    10930   23
chr1    10933   23
chr1    10936   1
chr1:10884
cd /home/shg047/oasis/biomark
less 7-P-3.PileOMeth.bedGraph_CpG.bedGraph
chr1    10857   10858   100     11      0
chr1    10860   10861   81      9	2
chr1    10866   10867   100     11      0
chr1    10884   10885   54      6	5
chr1    10886   10887   36      4	7
chr1    10902   10903   100     11      0
chr1    10907   10908   100     11      0
chr1    10913   10914   90      10      1
chr1    10923   10924   90      9	1
chr1    10928   10929   90      9	1
chr1    10930   10931   80      8	2

cd /home/shg047/oasis/monod/bam

samtools view -h -b 7-P-3.sorted.clipped.bam chr1:10854-10933 -o hg19.bsseq.chr1.10854-10933.bam
samtools tview -p chr1:10854 hg19.bsseq.chr1.10854-10933.bam ~/oasis/db/hg19/meth/bismark/hg19.fa
samtools mpileup -A -Q 1 --reference ~/oasis/db/hg19/meth/bismark/hg19.fa hg19.bsseq.chr1.10854-10933.bam
PileOMeth extract  -q 0 -p 1 /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa hg19.bsseq.chr1.10854-10933.bam  -o hg19.bsseq.chr1.10854-10933.bam.PileOMeth.bedGraph

PileOMeth extract --mergeContext --keepSingleton --keepDiscordant /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa hg19.bsseq.chr1.10854-10933.bam  -o hg19.bsseq.chr1.10854-10933.bam.PileOMeth.bedGraph



# 2016-2-10
Getopt::Std
  getopt  	
     getopt ('lw'); $opt_l and $opt_w to accept the value from terminal 
	 getopts ('abl:w:'); 

Getopt::Long
   GetOptions ('a|all' => \$all, 'l|length=i' => \$length,'w|width=i' => \$width);
     * for the option words, a double dash is required: ‘--length 24’ is acceptible
     * reference of $varabile should be taken as destination
     * You do not need to specified the option destination. If no destination is specified, GetOptions will define variables $opt_xxx where xxx is the name of the option, just like getopt and getopts. GetOptions will also accept a reference to a hash as its first argument and deliver the option values there, again just like getopt and getopts.
   GetOptions ('foo=i' => \@values);
     * Calling this program with arguments ‘-foo 1 -foo 2 -foo 3’ will result in @values having the value (1,2,3) provided it was initially empty.
     * Also, the option destination can be a reference to a hash. In this case, option values can have the form ‘key=value’. The value will be stored in the hash with the given key.


# 2016-2-10
PileOMeth extract  -q 10 -p 5 --minDepth 5 /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam  -o 7-P-3.bedGraph 
PileOMeth extract  -q 10 -p 5 --minDepth 5 -l biomark.list.bed /home/shg047/oasis/db/hg19/meth/bismark/hg19.fa /oasis/tscc/scratch/shg047/monod/bam/7-P-3.sorted.clipped.bam  -o  a

perl pipeline.pl biomark.list.bed /oasis/tscc/scratch/shg047/monod/bam/


bowtie2 -a -x ~/db/aligndb/hg19/bismark/Bisulfite_Genome/CT_conversion/BS_CT  ZNF154.fastq
bowtie2 -a -x ~/db/aligndb/hg19/bismark/Bisulfite_Genome/GA_conversion/BS_GA  ZNF154.fastq



BS_CT.rev.1.bt2
 perl -p -i -e 's/--paired-end/--single-end/g' *job

 
/oasis/tscc/scratch/ddiep/Working/Rerun_rrbs/BAMfiles/*


for i in `ls *methylFreq`
do
perl methylFreq2wig.pl $i 5 < $i > $i.wig 
gzip $i.wig
done


# 2016-2-5
@ZNF154-FP
GGTTTTTATTTTAGGTTTGA
+
aaaaaaaaaaaaaaaaaaaa
@ZNF154-RP
AAATCTATAAAAACTACATTACCTAAAATACTCTA
+
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
cd /home/shg047/work/biomarker/ZNF154
bismark --bowtie2 --non_directional --ambiguous --phred64-quals --fastq -L 10 -N 1 -s 0 --multicore 1 /home/shg047/db/aligndb/hg19/bismark ZNF154.
# 2016-2-5
bwa mem -O 0 -R '@RG    ID:ON1_S1	SM:ON1  PL:ILLUMINA     LB:On_S1' /home/zhl002/ZL_LTS33/RNA_seq/keiOTS_01262016/CAGOn-Cas_S1.fa 1-CAGOn-Cas_S1_L001_R1_001.fastq > CAGOnS1.fastq.sam
/home/kunzhang/softwares/samtools-latest/samtools view -b -t /home/zhl002/ZL_LTS33/RNA_seq/keiOTS_01262016/CAGOn-Cas_S1.fa.fai -h CAGOnS1.fastq.sam > CAGOnS1.fastq.bam
java -Xmx4g -jar /home/shg047/software/picard-tools-1.113/SortSam.jar TMP_DIR=./tmp/ INPUT=CAGOnS1.fastq.bam OUTPUT=CAGOnS1.sorted.bam QUIET=True SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
/home/kunzhang/softwares/samtools-latest/samtools index CAGOnS1.sorted.bam
java -Djava.io.tmpdir=./tmp -Xmx4g -jar /home/kunzhang/softwares/GenomeAnalysisTK-3.3/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/zhl002/ZL_LTS33/RNA_seq/keiOTS_01262016/CAGOn-Cas_S1.fa -I CAGOnS1.sorted.bam -glm INDEL -stand_call_
CAGOnS1.sh (END)
find -type f -iname '*.gz' -exec gunzip -t {}
/bsmap-2.74/bsmap -u -s 12 -v 0.04 -p 4 -a downsampled.SRR1232303_1.fastq -d hg19.fa -o downsampled.SRR1232303_1.fastq.sam
bsmap -u -s 12 -v 0.04 -p 6 -a ../fastq_trim/SRR1232304_1_trimmed.fq.gz -b ../fastq_trim/SRR1232304_2_trimmed.fq.gz -d  /home/shg047/oasis/Xliu2014/fastq -o SRR1232304_2_bam 
# 2016-2-4
cd /home/shg047/monod/hapinfo 
perl cgLD_Analysis_haploInfo_permuteRsq_v2_makeLD-plot.pl result.rsq chr10:101089204-101089305 < RRBS-7P30.hapInfo.txt
perl cgLD_Analysis_haploInfo_permuteRsq_v2_makeLD-plot.pl result.rsq chr10:101294771-101294802  < RRBS-7P30.hapInfo.txt
perl cgLD_Analysis_haploInfo_permuteRsq_v2_makeLD-plot.pl result.rsq chr10:103113834-103114495  < RRBS-7P30.hapInfo.txt
chr10:100995942-100996028

SRR1232310_1_val_1.fq.gz.temp.2_bismark_bt2_PE_report.txt
SRR1232313_1_val_1.fq.gz.temp.6_bismark_bt2_PE_report.txt
# 2016-1-27
Pregnancy.10.read1.fq.gz.job:samtools sort ../bam/Pregnancy.10.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.10.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.11.read1.fq.gz.job:samtools sort ../bam/Pregnancy.11.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.11.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.1.read1.fq.gz.job:samtools sort ../bam/Pregnancy.1.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.1.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.2.read1.fq.gz.job:samtools sort ../bam/Pregnancy.2.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.2.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.3.read1.fq.gz.job:samtools sort ../bam/Pregnancy.3.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.3.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.4.read1.fq.gz.job:samtools sort ../bam/Pregnancy.4.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.4.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.5.read1.fq.gz.job:samtools sort ../bam/Pregnancy.5.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.5.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
Pregnancy.6.read1.fq.gz.job:samtools sort ../bam/Pregnancy.6.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/Pregnancy.6.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
T21.3.read1.fq.gz.job:samtools sort ../bam/T21.3.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/T21.3.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
T21.4.read1.fq.gz.job:samtools sort ../bam/T21.4.read1_val_1.fq.gz_bismark_bt2_pe.bam -o ../bam/T21.4.read1_val_1.fq.gz_bismark_bt2_pe.sort.bam
# 2016-1-27

/home/kunzhang/CpgMIP/MONOD/Public_data/find_LMS.pl > whole_blood_LMS.txt
./extract_clusters.pl whole_blood_LMS.txt > whole_blood_LMS_clusters.txt

./report_tumor_HMH_regions.pl UCSD-004-05_sampleInfo.txt  | sort -k 3nr > UCSD-004-05_HMH_regions.txt
./report_tumor_HMH_regionsreport_tumor_HMH_regions.pl UCSD-004-07_sampleInfo.txt  | sort -k 3nr > UCSD-004-07_HMH_regions.txt
  
cat s_1_1_ILMN_Indx04.hapInfo.txt s_1_1_ILMN_Indx05.hapInfo.txt s_1_2_ILMN_Indx04.hapInfo.txt s_1_2_ILMN_Indx05.hapInfo.txt > PC-T-2.pooled.hapInfo.txt     
cat s_2_1_ILMN_Indx02.hapInfo.txt s_2_2_ILMN_Indx02.hapInfo.txt > PC-P-2.pooled.hapInfo.txt      
cat s_2_1_ILMN_Indx14.hapInfo.txt s_2_2_ILMN_Indx14.hapInfo.txt s_2_1_ILMN_Indx15.hapInfo.txt s_2_2_ILMN_Indx15.hapInfo.txt s_2_1_ILMN_Indx16.hapInfo.txt s_2_2_ILMN_Indx16.hapInfo.txt s_2_1_ILMN_Indx27.hapInfo.txt s_2_2_ILMN_Indx27.hapInfo.txt NP-RRBS-NC-P-1ng-p2-Jun20-A013.hapInfo.txt > NC-plasma.pooled.hapInfo.txt
 
../prepare_plotting_files.pl NC-plasma.pooled.hapInfo.txt PC-T-2.pooled.hapInfo.txt PC-P-2.pooled.hapInfo.txt chr6:152128536-152129155
../prepare_plotting_files.pl NC-plasma.pooled.hapInfo.txt PC-T-2.pooled.hapInfo.txt PC-P-2.pooled.hapInfo.txt chr5:176543913-176544076

# 2016-1-20

Cython: http://cython.org/
tar xzvf Cython-0.23.4.tar.gz
cd ./Cython-0.23.4/
sudo python setup.py install

wheel: https://pypi.python.org/pypi/wheel
tar xzvf wheel-0.26.0.tar.gz
cd ./wheel-0.26.0/
sudo python setup.py install

setuptools: https://packaging.python.org/en/latest/projects/#setuptools
tar xzvf setuptools-19.4.zip
cd ./setuptools-19.4
sudo python setup.py install

NumPy: http://www.scipy.org/scipylib/download.html
git clone git://github.com/numpy/numpy.git numpy
cd ./numpy/
sudo python setup.py install

MethylPurify: https://pypi.python.org/pypi/MethylPurify
tar xzvf MethylPurify-2.0-20141116.tar.gz
cd MethylPurify-2.0-20141116/
sudo python setup.py install --user
MethylPurify 


/bsmap-2.74/bsmap -u -s 12 -v 0.04 -p 4 -a downsampled.SRR1232303_1.fastq -d hg19.fa -o downsampled.SRR1232303_1.fastq.sam
python2.7 /home/ddiep/Downloads/MethylPurify-2.0-20140819/methylpurify/bin/MethylPurify -f downsampled.SRR1232303_1.fastq.bam 
      -g hg19.fa -i /home/ddiep/Downloads/MethylPurify-2.0/methylpurify/db/CGI_hg19_slop1000.bed -c 10 -s 50 -b 300
 
 
compare the difference between new and old perl script of GATK
old：/home/ajgore/AG_Ext12T/GATK_01022012/variantCallerBwaGATK-latest
new：/home/kunzhang/bin/variantCallerBwaGATK_05012015.pl


# 2016-1-19
bismark.zero.cov

ENCFF000LUN_trimmed

chr8:22249850-22249898

samtools tview -p chr8:22249850 /home/shg047/oasis/Haib/bam/ENCFF000LUN_trimmed.fq.gz_bismark_bt2.sort.bam /home/shg047/db/hg19.fa

qsub Lymphoma.run1.*job
qsub Pregnancy.14.*job
qsub T21.2.read1.*job

-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 07:17 Lymphoma.run1.read1_val_1.fq.gz.temp.4_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  453 Jan 19 07:19 Lymphoma.run1.read1_val_1.fq.gz.temp.5_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group 3.7K Jan 19 07:19 Lymphoma.run1.read1_val_1.fq.gz.temp.4_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group  453 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.6_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.5_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group  23K Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.6_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.1_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.3_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 07:20 Lymphoma.run1.read1_val_1.fq.gz.temp.2_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 09:20 Pregnancy.14.read1_val_1.fq.gz.temp.5_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 09:20 Pregnancy.14.read1_val_1.fq.gz.temp.3_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 09:21 Pregnancy.14.read1_val_1.fq.gz.temp.4_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  451 Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.6_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group 146K Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.3_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.6_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.4_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group  25K Jan 19 09:23 Pregnancy.14.read1_val_1.fq.gz.temp.5_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:18 T21.2.read1_val_1.fq.gz.temp.4_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:18 T21.2.read1_val_1.fq.gz.temp.3_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:18 T21.2.read1_val_1.fq.gz.temp.2_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:19 T21.2.read1_val_1.fq.gz.temp.1_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  437 Jan 19 11:19 T21.2.read1_val_1.fq.gz.temp.5_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group   62 Jan 19 11:19 T21.2.read1_val_1.fq.gz.temp.5_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group    0 Jan 19 11:19 T21.2.read1_val_1.fq.gz.temp.6_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  22K Jan 19 11:21 T21.2.read1_val_1.fq.gz.temp.2_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.2G Jan 19 17:39 T21.2.read1_val_1.fq.gz.temp.1_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.2G Jan 19 17:40 T21.2.read1_val_1.fq.gz.temp.4_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 2.3G Jan 19 17:40 T21.2.read1_val_1.fq.gz.temp.3_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 2.3G Jan 19 17:40 T21.2.read1_val_1.fq.gz.temp.6_bismark_bt2_pe.bam


# 2015-12-27
qsub Lymphoma.run1.read1.fq.gz.job
qsub Lymphoma.run2.read1.fq.gz.job
qsub Lymphoma.run3.read1.fq.gz.job
qsub Pregnancy.14.read1.fq.gz.job
qsub Pregnancy.9.run2.read1.fq.gz.job
qsub T21.2.read1.fq.gz.job

rm Pregnancy.13.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.12.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.15.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.7.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.8.run2.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.9.run2.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.8.run1.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.9.run1.read1_trimmed.fq.gz_bismark_bt2_*
rm Pregnancy.14.read1_trimmed.fq.gz_bismark_bt2_*
rm T21.2.read1_trimmed.fq.gz_bismark_bt2_*
rm T21.5.read1_trimmed.fq.gz_bismark_bt2_*
rm T21.1.read1_trimmed.fq.gz_bismark_bt2_*


ls Lymphoma.run3* | grep -v val 

ls *pe 

qsub HL_05.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_06.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_07.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_08.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_09.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_14.sorted.clipped.fastq.gz.bismark.pbs
qsub HL_15.sorted.clipped.fastq.gz.bismark.pbs

ls *LTP* | grep -v val

rm LTP2.read1_trimmed.fq.gz_bismark_bt2_pe.bam
rm LTP2.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
rm LTP3.read1_trimmed.fq.gz_bismark_bt2_pe.bam
rm LTP3.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
rm LTP4.read1_trimmed.fq.gz_bismark_bt2_pe.bam
rm LTP4.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt

ls -larth *Pregnancy* | grep -v val

-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 00:30 Pregnancy.5.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  21G Jan 15 00:30 Pregnancy.5.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 02:03 Pregnancy.4.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  23G Jan 15 02:03 Pregnancy.4.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group  17G Jan 15 02:45 Pregnancy.10.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 02:45 Pregnancy.10.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 04:04 Pregnancy.6.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  22G Jan 15 04:04 Pregnancy.6.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 04:36 Pregnancy.2.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  26G Jan 15 04:36 Pregnancy.2.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 07:23 Pregnancy.3.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  30G Jan 15 07:23 Pregnancy.3.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 09:03 Pregnancy.11.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  24G Jan 15 09:03 Pregnancy.11.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 11:15 Pregnancy.1.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  31G Jan 15 11:15 Pregnancy.1.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:42 Pregnancy.13.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  707 Jan 15 13:42 Pregnancy.13.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:42 Pregnancy.12.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  707 Jan 15 13:42 Pregnancy.12.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:44 Pregnancy.15.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  707 Jan 15 13:44 Pregnancy.15.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:44 Pregnancy.7.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  706 Jan 15 13:44 Pregnancy.7.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.9K Jan 15 13:44 Pregnancy.8.run2.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  710 Jan 15 13:44 Pregnancy.8.run2.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.9K Jan 15 13:44 Pregnancy.9.run2.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  710 Jan 15 13:44 Pregnancy.9.run2.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.9K Jan 15 13:44 Pregnancy.8.run1.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  710 Jan 15 13:44 Pregnancy.8.run1.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.9K Jan 15 13:44 Pregnancy.9.run1.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  710 Jan 15 13:44 Pregnancy.9.run1.read1_trimmed.fq.gz_bismark_bt2_pe.bam
-rw-r--r-- 1 shg047 k4zhang-group 1.8K Jan 15 13:45 Pregnancy.14.read1_trimmed.fq.gz_bismark_bt2_PE_report.txt
-rw-r--r-- 1 shg047 k4zhang-group  706 Jan 15 13:45 Pregnancy.14.read1_trimmed.fq.gz_bismark_bt2_pe.bam



# back-up data
# WGBS MHL matrix(update by Shicheng) were save in (UTH):
/home/sguo/monod/hap/wgbs/All_chromosomes_combined/wgbs.mhl.txt
# RRBS(Batch 1 and 2) MHL matrix(update by Shicheng) were save in (UTH):
/home/sguo/monod/rrbs_kun/rrbs.kun.mhl.txt
/home/sguo/monod/rrbs_encode
cd /oasis/tscc/scratch/k4zhang/MONOD
cp /oasis/tscc/scratch/k4zhang/MONOD/Ecker_Tissue_WGBS/*pl ./
cp /oasis/tscc/scratch/k4zhang/MONOD/tumor_WGBS/*pl ./
TSCC_whole_blood_WGBS_sampleInfo_WGBS-pooled-mld-blocks.txt
cp /oasis/tscc/scratch/k4zhang/MONOD/whole_blood_WGBS/batch_bam2hapInfo2.pl ./
# 2015-12-20
wc -l /home/shg047/db/hg19/hg19_refGene.bed
head /home/shg047/db/hg19/hg19_refGene.Fantom.bed
cd ~/bioin/annotation/
perl /home/shg047/bioin/bin/trim.pl
perl /home/shg047/bioin/bin/select.pl
perl /home/shg047/bioin/bin/unique.pl
less /home/shg047/bioin/bin/unique.pl
# build new hg19_reference, replace enhancer with Fantom Enhancer
cd ~/db/hg19
grep -v Enhancer hg19_refGene.bed > hg19_refGene.2.bed
awk '{print $1,$2,$3,$4,$5,$6,"Enhancer",$8,$9}' OFS="\t" Enhancers.Fantom.hg19.bed > Enhancer.2.Fantom.hg19.bed9
cat Enhancer.2.Fantom.hg19.bed9 hg19_refGene.2.bed > hg19_refGene.Fantom.bed
rm Enhancer.2.Fantom.hg19.bed9
rm hg19_refGene.2.bed

cd ~/monod/methyblock
# loop source files and get all the overlap(A,B) and non-overlap(A)regions
java -jar ~/bin/M2G.jar -s ~/db/hg19/hg19_refGene.Fantom.bed -l WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed -o rlt
# select overlap regions with source by 9th column
perl /home/shg047/bioin/bin/trim.pl rlt 9
perl /home/shg047/bioin/bin/select.pl rlt.true 6
perl /home/shg047/bioin/bin/unique.pl rlt.true

# 2015-12-18
cat GSM1051152_colon.CpG.calls.txt.bedgraph | awk 'NR<2{print $0;next}{print $0| "sort -k1,1 -k2,2n"}' >  GSM1051152_colon.CpG.calls.txt.sort.bedgraph &
samtools view -h -o human.aligned.cleaned.sam song.aligned.cleaned.bam
samtools view -h -o mouse.aligned.cleaned.sam mouse.aligned.cleaned.sam

# 2015-12-17
grep -v enhancer enhancer.12.encode.cell.type.txt | grep -v file | grep -v Chromosome > enhancer.12.encode.cell.types.txt
perl -lane "next if ! /^chr/; print" enhancer.12.encode.cell.types.txt > enhancer.12.encode.cell.types.true.txt
dos2unix enhancer.12.encode.cell.types.true.txt 
perl -lane "s/\s/\t/g" enhancer.12.encode.cell.types.true.txt 
awk '{print $1, $2-1000, $2+1000,$3}' OFS="\t" enhancer.12.encode.cell.types.true.txt > enhancer.12.encode.cell.types.hg18.txt 
perl -lane "next if @F[0]<0, print" enhancer.12.encode.cell.types.hg18.txt > enhancer.12.encode.cell.types.hg18.nonneg.txt 
wc -l enhancer.12.encode.cell.types.hg18.txt
wc -l enhancer.12.encode.cell.types.hg18.nonneg.txt
bedtools sort -i enhancer.12.encode.cell.types.hg18.nonneg.txt > enhancer.12.encode.cell.types.hg18.sort.txt 
perl -lane "print if @F[-1]>0.95" enhancer.12.encode.cell.types.hg18.sort.txt > enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt
bedtools merge -i enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt > enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt
wc -l enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt
wc -l enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt
lifeOver

# 2015-12-10
/home/shg047/monod/methyblock/wgbs
wc -l /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/All_WGBS_pooled/mappable_bins_hapInfo/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
wc -l wgbs.bed
bedtools intersect -wao -a wgbs.bed -b ~/db/hg19/CpGI.hg19.bed | perl -lane '{print if @F[8]==0}' | sort -u | wc -l 
cp /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/All_WGBS_pooled/mappable_bins_hapInfo/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed
bedtools intersect -v -a WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed -b ~/db/hg19/CpGI.hg19.bed | wc -l

WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed

cd 
/home/shg047/software/broadinstitute-picard-2a49ee2/src/java/picard/sam
/home/songchen/RNASeq/151211/fastq
scp /media/LTS_33T/SeqStore2016/151211_MiSeq/Fastq
git agordon/libgtextutils
/home/shg047/software/bin/fastx_trimmer

scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/STAR_Index_Gencode_75bp  shg047@genome-miner.ucsd.edu:/home/shg047/db/song
scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/hg38.fa  shg047@genome-miner.ucsd.edu:/home/shg047/db/song
scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/gencode.v23.annotation.gtf  shg047@genome-miner.ucsd.edu:/home/shg047/db/song

scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/STAR_Index_Mouse_Gencode  shg047@genome-miner.ucsd.edu:/home/shg047/db/song
scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/GRCm38.fa  shg047@genome-miner.ucsd.edu:/home/shg047/db/song
scp -r /oasis/tscc/scratch/yaw004/Genome_Indexes/gencode.vM7.annotation.gtf  shg047@genome-miner.ucsd.edu:/home/shg047/db/song

mouseGenomeIndex=/oasis/tscc/scratch/yaw004/Genome_Indexes/STAR_Index_Mouse_Gencode
mouseReferenceFasta=
mouseReferenceGTF=


cp /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/All_WGBS_pooled/mappable_bins_hapInfo/no_tumor/WGBS_pooled_mappable_bins.mld_blocks_r2-0.5.bed  /home/shg047/monod/methyblock

source("/home/shg047/monod/phase3/bedEnrichmentTest.R")



/home/shg047/db/tmp/liftOver

perl -lane "next if ! /^chr/; print" enhancer.12.encode.cell.types.txt > enhancer.12.encode.cell.types.true.txt
dos2unix enhancer.12.encode.cell.types.true.txt 
perl -lane "s/\s/\t/g" enhancer.12.encode.cell.types.true.txt 
awk '{print $1, $2-1000, $2+1000,$3}' OFS="\t" enhancer.12.encode.cell.types.true.txt > enhancer.12.encode.cell.types.hg18.txt 
perl -lane "next if @F[0]<0, print" enhancer.12.encode.cell.types.hg18.txt > enhancer.12.encode.cell.types.hg18.nonneg.txt 
wc -l enhancer.12.encode.cell.types.hg18.txt
wc -l enhancer.12.encode.cell.types.hg18.nonneg.txt
bedtools sort -i enhancer.12.encode.cell.types.hg18.nonneg.txt > enhancer.12.encode.cell.types.hg18.sort.txt 
liftOver enhancer.12.encode.cell.types.hg18.sort.txt  ~/bin/hg18ToHg19.over.chain.gz enhancer.12ct.all.hg19.txt  tmp
cp enhancer.12ct.all.hg19.txt ~/db/hg19/enhancer.12ct.encode.all.hg19.txt

perl -lane "print if @F[-1]>0.95" enhancer.12.encode.cell.types.hg18.sort.txt > enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt
bedtools merge -i enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt > enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt
wc -l enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.txt
wc -l enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt
liftOver /home/shg047/monod/RRBS_Plasma_Batch2/enhancer.12.encode.cell.types.hg18.sort.cutoff0.95.merge.txt ~/bin/hg18ToHg19.over.chain.gz enhancer.12ct.encode.hg19.bed tmp

bedtools window -u -w 0 -b enhancer.12ct.all.hg19.txt -a ~/db/hg19/Enhancers.Fantom.hg19.bed | wc -l 



# 2015-12-10
bedtools shuffle -i /home/shg047/monod/methyblock/wgbs.bed -incl ~/db/hg19/hg19_refGene.bed -g ~/db/hg19/hg19.chrom.sizes 
bedtools intersect -wo -a /home/shg047/monod/methyblock/wgbs.bed -b ~/db/hg19/hg19_refGene.bed 

perl bedAnnoEnrich.pl -i /home/shg047/monod/methyblock/wgbs.bed -r ~/db/hg19/hg19_refGene.bed -g ~/db/hg19/hg19.chrom.sizes -n 50

# compute mhl from haploinfo files
perl /oasis/tscc/scratch/shg047/RRBS/src/haplo2hml.pl /oasis/tscc/scratch/shg047/RRBS/batch2/hapinfo > ../rrbs.batch2.mhl.txt

find ./ -type f -name *pl -exec chmod 700


# copy mhl file from tscc to genome-miner
scp rrbs.batch2.mhl.txt shg047@genome-miner.ucsd.edu:/home/shg047/monod/RRBS_Plasma_Batch2

# R 




cd /oasis/tscc/scratch/shg047/RRBS/
cd /home/shg047/db/hg19
cd /oasis/tscc.old/scratch/k4zhang/MONOD
cd /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles
cd /home/ddiep/Plasma_RRBS/BAMfiles


cd /oasis/tscc/scratch/ddiep/Plasma_RRBS_151208/BAMfiles
scp shicheng@meangenemachine.dynamic.ucsd.edu:/home/ddiep/Plasma_RRBS/BAMfiles/RRBS-6P24.sorted.clipped.bam ./
scp shicheng@meangenemachine.dynamic.ucsd.edu:/home/ddiep/Plasma_RRBS/BAMfiles/RRBS-6P24.sorted.clipped.bam.bai ./



awk '{print $1}' combined_list_files | sort -u 
 
less /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed

cp /oasis/tscc.old/scratch/k4zhang/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed ./

scp shicheng@meangenemachine.dynamic.ucsd.edu:



/home/shg047/db/hg19/HsGenome19.CpG.positions.txt

scp shg047@genome-miner.ucsd.edu:/home/shg047/db/hg19/HsGenome19.CpG.positions.txt ./



/home/shg047/db/hg19/genomeDistribution
wget -qO- ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gff3.gz \
    | gunzip --stdout - \
    | gff2bed - \
    > annotations.bed

# 2015-12-3

screen -ls | grep Detached | cut -d. -f1 | awk '{print $1}' | xargs kill
screen -ls | grep pts | cut -d. -f1 | awk '{print $1}' | xargs kill


# 2015-12-1


for i in `ls  methylC*.txt.pairwiseR2`
do
cp pairwise.R $i.pairwise.R
R CMD BATCH "--args $i" $i.pairwise.R & 
done

for i in `ls  N37*C*.txt.pairwiseR2`
do
cp pairwise.R $i.pairwise.R
R CMD BATCH "--args $i" $i.pairwise.R & 
done


R CMD BATCH "--args WB_new-born.all_chrs.hapInfo.txt.pairwiseR2" WB_new-born.all_chrs.hapInfo.txt.pairwiseR2.pairwise.R & 





R CMD BATCH "--args N37-Liver.all_chrs.hapInfo.txt.pairwiseR2" $i.pairwise.R & 
R CMD BATCH "--args $i" $i.pairwise.R & 
R CMD BATCH "--args $i" $i.pairwise.R & 

for i in `ls *hapInfo.txt.pairwiseR2`
do
cp pairwise.R $i.pairwise.R
R CMD BATCH "--args $i" $i.pairwise.R &
done



for i in `ls /home/shg047/db/hg19/encode*bed`
do
echo $i
bedtools intersect -wb -a en.bed -b $i >> en.encode.share.bed
done



for i in `ls /home/shg047/db/hg19/encode*bed`
do
echo $i
bedtools intersect -wb -a en.bed -b $i >> en.encode.share.bed
done

for i in `ls /home/shg047/db/hg19/encode*bed`
do
echo $i
bedtools intersect -wb -a ec.bed -b $i >> ec.encode.share.bed
done

for i in `ls /home/shg047/db/hg19/encode*bed`
do
echo $i
bedtools intersect -wb -a me.bed -b $i >> me.encode.share.bed
done


awk '{print $7}' en.encode.share.bed | sort -u > en.encode.tf.txt
awk '{print $7}' ec.encode.share.bed | sort -u > ec.encode.tf.txt
awk '{print $7}' me.encode.share.bed | sort -u > me.encode.tf.txt





cd /oasis/tscc/scratch/shg047/N37
samtools merge all.bam `find /basedir/ -name "*myfiles*.bam"`

samtools view -H  Indx01.chr6.rmdup.bam > /oasis/tscc/scratch/shg047/N37/merge.sam.header
 cd /home/k4zhang/my_oasis_tscc/MONOD/N37_WGBS/BAMfiles
samtools merge -h /oasis/tscc/scratch/shg047/N37/merge.sam.header /oasis/tscc/scratch/shg047/N37/Indx01.bam Indx01.chr1.rmdup.bam Indx01.chr10.rmdup.bam 

samtools view -H Indx01.chr1.rmdup.bam
samtools view -H Indx01.chr2.rmdup.bam


bam2haploinfo.pl
64 to phred

WGBS_pooled_mappable_bins.chr1.mld_blocks_r2-0.5.bed

1, The refGene.txt file is a database file, and consequently is based on the internal representation.
2, Our internal database representations of coordinates always have a zero-based start and a one-based end. We add 1 to the start before displaying coordinates in the Genome Browser. 
3, If you submit data to the browser in position format (chr#:##-##), the browser assumes this information is 1-based. 
4, If you submit data in any other format (BED (chr# ## ##) or otherwise), the browser will assume it is 0-based
5, Similarly, any data returned by the browser in position format is 1-based, while data returned in BED format is 0-based.
6, Some file formats are 1-based (GFF, SAM, VCF) and others are 0-based (BED, BAM)


cd /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/
samtools view -q 20 methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam chr1:10003741-10003773 > methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.sam

/home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam
/home/shg047/monod/haplo/mergedBam2hapInfo_WGBS_25Jun15.pl /home/k4zhang/my_oasis_tscc/MONOD/All_WGBS_pooled/WGBS_pooled_mappable_bins.chr1.mld_blocks_r2-0.5.bed  /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam > methylC-seq_h1-npc_r1.chr1.10003741.10003773.hapInfo.txt

/home/kunzhang/CpgMIP/MONOD/Data/141216_HiSeqRapidRunSeqCap/mappable_bin_hapInfo
less 141216_SeqCap_tumor_RD10_80up.chr21.hapInfo.txt

 head /home/k4zhang/human_g1k_v37/HsGenome19.CpG.positions.txt
 

grep /home/k4zhang/human_g1k_v37/HsGenome19.CpG.positions.txt chr1\s


# 2015-11-30
/home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/All_chromosomes_combined
cd /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/
samtools view -b /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/methylC-seq_h1-npc_r1.chr1.rmdup.bam chr1:10003741-10003773 > methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam
samtools index methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam

cd 
/home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles
mv methylC-seq_h1-npc_r1.chr1.10003741.10003773.rmdup.bam /home/shg047/monod/



# 2015-11-20

1, download the data
2, checkphred
3, trim_galore
4, compress to gz
5, fastqc
/home/k4zhang/softwares/STAR_2.3.0e.Linux_x86_64/STAR --readFilesCommand zcat --runThreadN 8  --outSAMstrandField intronMotif   --genomeDir /home/k4zhang/my_oasis_tscc/StarIndex/hg19 --readFilesIn Sample_AL-PGP1ipstoCM-102315-3_S3/AL-PGP
# db for bismarker alignment

cd ~/db/aligndb/hg19/bismark
scp hg047@genome-miner.ucsd.edu:/home/shg047/db/hg19/hg19.fa /home/shg047/db/hg19/meth/bismark/Bisulfite_Genome


bismark --bowtie2 --phred64-quals --fastq  --non_directional
bismark -q --phred33-quals --non_directional --bowtie2 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq

cd /oasis/tscc/scratch/shg047/DennisLo2015/fastq

bismark --bowtie2 --phred64-quals --fastq -L 30 -N 1 /home/shg047/db/hg19/meth/bismark -1 Lymphoma.run2.read1.fq.gz -2 Lymphoma.run2.read2.fq.gz -o ../bam




# fastq dataset
cd /oasis/tscc/scratch/shg047/DennisLo2015/fastq

cd /media/TmpStore1/DennisLo2015/


# FastQC for large numbers of samples
http://www.cureffi.org/2012/08/27/fastqc-for-large-numbers-of-samples/?utm_source=delicious+via+twitterfeed&utm_medium=twitter

# 2015-11-17
backtick · bash cmd ·

set filename but leave passwd empty and then 
ssh-keygen
ssh-copy-id shg047@genome-miner.ucsd.edu
Install perl module

wget https://s3.amazonaws.com/deqiangsun/software/moabs/moabs-v1.3.2.src.x86_64_Linux.data.tar.gz
tar xzvf moabs-v1.3.2.src.x86_64_Linux.data.tar.gz
cd moabs-v1.3.2.src.x86_64_Linux.data
cpan App::cpanminus
cpanm Config::Simple
make
pwd
export PATH=/home/shg047/software/moabs-v1.3.2.src.x86_64_Linux.data/bin/:$PATH

moabs --cf mytestrun.cfg
cd data/
moabs -i wt_r1.fq -i wt_r2.fq -i ko_r1.fq -i ko_r2.fq
mcomp -r wt_r1.bam.G.bed,wt_r2.bam.G.bed -r ko_r1.bam.G.bed,ko_r2.bam.G.bed -m wildtype -m knockout -c comp.wiVar.txt --withVariance 1 -p 4 



cd /oasis/tscc/scratch/shg047/DennisLo2015
Ricky Martin

#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;
my @file=glob("*.phred64.fq");
foreach my $file(@file
my ($sample,undef)=split /\.phred64.fq/,$file;
rename $file "$sample.fq";
system("gzip $sample.fq");
print "$file\t$phred\n";
}



for i in `*fq`
do
reformat.sh in=$i out=$i.phred64.fq qin=33 qout=64
done



# 2015-11-18
cp /home/kunzhang/bin/trimFastq.pl ~/bin/
vim  ~/bin/trimFastq.pl
vim  /home/kunzhang/bin/trimFastq.pl
cd /home/zhl002/ZL_LTS33/RNA_seq/111215_PGP1iPstoCM/data
perl ~/bin/trimFastq.pl AL-PGP1ipstoCM-102315-10_S10_L001_R1_001.fastq
perl /home/kunzhang/bin/trimFastq.pl < AL-PGP1ipstoCM-102315-10_S10_L001_R1_001.fastq 0 20


ssh-keygen -t rsa
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDQig1C3Oz7sn8Z7dgYWjJL8oyZH71ul+EqmggMC6CWWRZNmXSx73kovKIDMUhrXTCiAXcQ+pHVr9N7PKewZ0Zd2IMfYWcr4tFC5orfOHL40pXyecm1tMwhFbv6U54iV3aPK5PGgGO8zy/n0Ckbj+gKWD95ySQPeptM4lFJ4dqp3z/TudIKAZJqOWj6VzG6MgUZTNpCXO76Lr8A4NXWmMumbcdK0LylSTBlCBjLfR61znP9sZpXEOZwv/OnNmr0pm7D4IM3ny+6sNXBtnexIVn3C15vRgA3OvyVHCZ89EsOLwzMxRcUU3Xp4KQVPNRmbQJ3Utfy9/uh6TE2EgFEZ+gf shg047@genome-miner2
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCtMUHApnb9L6GvhSgaQpbogwioZXjTpi3gg3yJg4YTrnvnfSyVoRz1ZBHkKYVW/HzGDOfLejVe2Pzbkyi06LoO05TiSTm6/ytV47jOMTMnuaJj/pDXMiGashYfy/Gm8TguYGFQeGtj5f+3V6hIGnIHcyF1YFf8qClDwgRQDv4DilM9FpTGU+HcuOIid9NuWcPEv96vwm4Z3mxjXLqAZS0X6aygU4ge6yGPl0zpl1aWX0lzcl0Q9c7mUZ8Z5elgUimX/Ogl+HhJapnhSNWudAV2UI7QXRUphcXdVjaVFRVAKSX+tC6bQ4vU2pjLOxHanadayM3mbpz+IIYaJFPw1B/n shg047@genome-miner2

vim ~/.ssh/id_rsa.pub

scp shg047@genome-miner.ucsd.edu:/media/TmpStore1/DennisLo2015/HOT227.fq.gz  /oasis/tscc/scratch/shg047/DennisLo2015/fastq/
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQC2G7U7lO5uw0a7ft9h1FnmU8ITeP/69A49QpHX+OJay+EE+yL6PGovOG4qbq/BnUXSh78u9FLnZxP/kb2ThAdwVz2Ihent2EhpzL4ka2tUirDe4mePUlqMCnyw9O+huwIkiUKUp+C/Wo7IaFs/bltmLpnz3cANIrrXvZr7r2+6nDtai2eZpCYDI1l2R6+cEGuZBA5fwJuXkX47tYHcD1g/ptFrKJzH8WJ9EurRzM8S8F7xeRJhAl3o57+NXzrGcrepEUHX7xfs+ofM3IHPjgVxpgA9875SjsSF4wvwh8XRzqo4FmF0LAv05LLlSSR9JZ4vwsNrM2KXwRw7CReyb+mr shg047@tscc-login1.sdsc.edu

scp shg047@genome-miner.ucsd.edu:/media/TmpStore1/DennisLo2015/HOT227.fq.gz  /oasis/tscc/scratch/shg047/DennisLo2015/fastq/

scp /media/TmpStore1/DennisLo2015/HOT227.fq.gz shg047@tscc-login.sdsc.edu:/oasis/tscc/scratch/shg047/DennisLo2015/fastq/
scp *.fq.gz shg047@tscc-login.sdsc.edu:/oasis/tscc/scratch/shg047/DennisLo2015/fastq/




# 2015-11-17
for i in `ls *fq.gz` 
do 
perl ~/bin/checkphred.pl $i
done

#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;
my @file=glob("*fq.gz");
foreach my $file(@file){
my $phred=`system(perl ~/bin/checkphred $file)`
print "$file\t$phred\n";
}


#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;
my $file=shift @ARGV;
open F, $file;





#!/usr/bin/bash
for i in `ls *fq.gz`
do
j=`perl ~/bin/checkphred.pl $i`
echo "$i $j"
done



grep ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCG AL-LAMPCRGFP-111015-N2indx16_S40_L001_R1_001.fastq | wc -l
grep CAAGTTCAGCGTGTCCG AL-LAMPCRGFP-111015-N2indx16_S40_L001_R1_001.fastq | wc -l


grep ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCG AL-LAMPCRGFP-111015-N2indx16_S40_L001_R1_001.fastq | cut -c 1-60 > keep.fa

perl ~/bin/fa2fastq.pl keep.fa > keep.fastq




# 2015-11-16
@HISEQ:565:C6NCYACXX:5:1101:1220:1868 1:N:0:
NGTTAAATGTTAGGTGAATTTTAGTTTTATTGTTGAGTGAGATTTTTTGCGTTTAGTGTTGTTTTTATATTAGTT
+
#1=DDD>DF?FFDFFFFFIIFIIIIICACFHHHHBDDEFGF>GFIFIFE907BFFFFG=FFIIIIFFEEFEE@=;
@HISEQ:565:C6NCYACXX:5:1101:1220:1868 2:N:0:
NGTTAAATGTTAGGTGAATTTTAGTTTTATTGTTGAGTGAGATTTTTTGCGTTTAGTGTTGTTTTTATATTAGTT
+
#1=DDD>DF?FFDFFFFFIIFIIIIICACFHHHHBDDEFGF>GFIFIFE907BFFFFG=FFIIIIFFEEFEE@=;

killall screen
screen -wipe
pkill screen




prepare fa and fastq for illumina adapter and primers
cd ~/software/FastQC/Configuration
perl txt2fastq.pl adapter_list.txt > adapter_primer.fastq
perl txt2fastq.pl contaminant_list.txt >> adapter_primer.fastq
perl txt2fasta.pl adapter_list.txt > adapter_primer.fa
perl txt2fasta.pl contaminant_list.txt >> adapter_primer.fa


cat read.list |awk '{print "trim_galore -q 20 --phred33 -a AGATCGGAAGAGC "$1}' >  Run.trimgalore.single.sh


# 2015-11-05
TopHat remove the part of read names which distinguish read1 and read2, thus, when using paired-end data as single-end with TopHat.

module load cuda
module load biotools

grep AGATCGGAAGAGC

gunzip -c T21.5.read1.fq.gz | head -n 400000 > T21.5.read1.fq.head
gunzip -c T21.5.read2.fq.gz | head -n 400000 > T21.5.read2.fq.head 

trim_galore --paired -a GATCGGAAGAGCA -a2 GCTCTTCCGATCT --retain_unpaired  --trim1  T21.5.read1.fq.head T21.5.read2.fq.head  

grep AGATCGGAAGAGC T21.5.read1.fq.head_unpaired_1.fq
grep GCTCTTCCGATCT T21.5.read2.fq.head_val_2.fq

grep AGATCGGAAGAGC T21.5.read1.fq.head
grep GCTCTTCCGATCT T21.5.read1.fq.head 
grep AGATCGGAAGAGC T21.5.read2.fq.head
grep GCTCTTCCGATCT T21.5.read2.fq.head 


gunzip -c T21.5.read1.fq.gz | grep AGATCGGAAGAGC | head
gunzip -c T21.5.read2.fq.gz | grep GCTCTTCCGATCT 

gunzip -c T21.5.read1.fq.gz | grep AGATCGGAAGAGC
gunzip -c T21.5.read2.fq.gz | grep TCTAGCCTTCTCG


trim_galore support gzip compressed FastQ files
bismark also support gzip compressed FastQ files (-input files can be uncompressed or gzip-compressed (ending in .gz))

gzip vs. bzip2
1, gzip is faster than bzip2 and the compress ratio is about 35% to 40%
2, bzip2 is slower than gzip however the compress ratio is about 34% to 35%
3, gzip is more conveniment in the usage since majority software can take gzip as the input files

# 2015-11-04


# 2015-11-03
cd /home/zhl002/ZL_LTS33/mouse_WGBS/PCR_validation
bismark -q --phred33-quals --non_directional --bowtie2 -n 1 -l 5 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq
bismark2report ./
deduplicate_bismark
coverage2cytosine
bismark2bedGraph
bismark_methylation_extractor


cd /media/TmpStore1/DennisLo2015
try to analysis the RRBS data from encode project and find the difference between encode project and our own data
1, coverage region
2, coverage

# encode project
cd /media/LTS_60T/shg047/encode
bismark -q --phred64-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark wgEncodeHaibMethylRrbsAg09309UwRawDataRep1.fastq &
bismark -q --phred64-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark wgEncodeHaibMethylRrbsAg09309UwRawDataRep2.fastq &
bismark -q --phred64-quals --bowtie2 -n 1 -l 10 /home/shg047/db/aligndb/hg19/bismark a.fastq 

grep TGTGAGATAGGTTTTAGGTGAAG *.fastq
grep GTTGTAAGAGTTAATGGTTTATTATAAT *.fastq
grep TTTGTTTTTTATTAAGTTAGATGTGT  *.fastq
grep GGTGTGAGTTTGTTTTTTTTAGTTTTTA  *.fastq
grep TTAGGAGGGTAAGGAATATTTTAGG *.fastq
grep TGTGAGTTGAAGTAGGAAGGTTTTT *.fastq
grep TTTGTTTGTTGTTTTTGGAGGATA *.fastq
grep GTTGATGTATGTTTTATTAGTAAATAA *.fastq
grep AGAATGTAGGTTTTTGTTTTGTAGA *.fastq
grep ACCCTAAATTATAATCTCAAACACC  *.fastq
grep AATATCTCTAACTCCTTAACCCACC  *.fastq
grep TCCAAATAAAACAATATTTAAAATCTCC   *.fastq
grep AAACAATTTACTTCAACACAACCAA   *.fastq

GCGAGAAGGCTAGTTTGTTTGTTGTTTTTGGAGGATAGAAGAGCGGTTCAGCAGGAATGCCGAGCAACCCATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAATTAATAACACAAATCAGAAACCCAAGTGACTGACAGAGAAACCC

/home/shg047/db/aligndb/mm9


cd /home/zhl002/ZL_LTS33/mouse_WGBS/PCR_validation
perl checkscore.pl  ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq ZL-MeValidate-enPCR-14-Oct21_S5_L001_R1_001.fastq ZL-MeValidate-enPCR-12-Oct21_S3_L001_R1_001.fastq

grep AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq
grep AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT *fastq

TGGCTATTTTAGCTGTACCTCAGTATGCCTTAGGGTCATCTTGGGTAATGGGTCTAGAGGACACTTTTAGTGCTAAAACCTGAAAGGGATCCAGTGGTGCCAATTTAGCAACCATGCAGGACAGTAGCCCTGCGGGGAATGACTTCT


grep GTAGGGTGGGAAGTGGTTTTTAT *fastq
grep GGTGTGAGTTTGTTTTTTTTAGTTTTTA *fastq
grep TTTGTTTGTTGTTTTTGGAGGATA *fastq
grep 

CGCAGACATTTCTTCCTGCTTCTCTGCTGTTCTCCGTTTTGTTGTCCTCTTCTTAGGCAGCAC
GCGAGAAGGCTAGTTTGTTTGTTGTTTTTGGAGGATAGAAGAGCGGTTCAGCAGGAATGCCGAGCAACCCATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAATTAATAACACAAATCAGAAACCCAAGTGACTGACAGAGAAACCC

TTGTTTTTAGTATAATTTTTTTATAGTGTGGTTTGGGTGGGGTTTTTAGAAGGTTTTTTTTGTTTTTAGAGTAAGGTTGATTTTTTTTAAGGTATTGTTTGTTTTTGGTTGTTTTTTATTGTTTTGTTTTTTATTAAGTTAGATGTGTGTTGTAGATATTTTTTTTTGTTTTTTTGTTGTTTTTTGTTTTGTTGTTTTTTTTTTAGGTAGTATGTAGGTTTTATTTGAATTGTTTTAGTATTGTTATGTTTTGGAGATTTTAAATATTGTTTTATTTGGAGTTTTAAATGGTTAAATTTTATGTGGTTTGAGGGTTGGTAGAAGTTATGTTTTTTTTTGGTTTTTTTTATTTTTTTTATGGTT
AATGAGTAATTTTAAGATAGTTTTTTAAAATTAAGATTATTTAGTGGTGGTTTTGTTTTGTTGTTTATTTTATGTTTGTGGGTGTTTTTTTTGTATGTATATATTTTTATTATGTATTTGTTTGTTGTTTTTGGAGGATAGAAGAGGGTATTGGAATTTTTTGGATTGGAGTTGTAGATGGTTGTGAGTTATTATATGGGTGTTTGGAGATGAGTTTAGGTTATTTTGAAGAGTAGTTATTGTATTTAATTGTTGAGTTATTTTTTGAGTTTTGGGTGTTTTGAGATAGGTAAGGGTATAGATGAAGATAATTTTTAAAGGTTTTTTGAGTTTTTTGGATTGTTATATGTTTTTT


fastqc -f fastq -t 3 T21.5.read2.fq.gz  # fastqc detect the file format automatically(fastq,bam,sam)
fastqc -t 6 T21.5.read2.fq.gz


trim_galore sequence_file.fastq
trim_galore --phred33 --fastqc --illumina --non_directional --rrbs --paired --dont_gzip -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq
trim_galore --phred33 --fastqc --non_directional --rrbs --dont_gzip ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq
trim_galore --phred33 --fastqc --non_directional --rrbs --dont_gzip -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq

bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq

samtools tview -p chr17:76728261 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
samtools tview -p chr8:12808189 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
samtools tview -p chr9:92831023 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
samtools tview -p chr5:27124338 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
samtools tview -p chr3:7632484 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa





cd /home/zhl002/ZL_LTS33/mouse_WGBS/PCR_validation
bismark -q --phred33-quals --non_directional --bowtie2 -n 1 -l 5 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq -o ../


samtools tview -p chr9:52407014 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm10/mm10.fa

samtools tview -p chr11:13303365 ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001_trimmed.fq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm10/mm10.fa


samtools tview -p chrX:115198611 ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq_bismark_bt2.sort.bam /home/shg047/db/aligndb/mm9/mm9.fa
chr11:13303365
chrX:115198611 
AGATCGGAAGAGC


fastq-mcf ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq
fastq-stats ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq

fastq-clipper
fastq-join
fastq-multx


trim_galore ZL-MeValidate-enPCR-19-Oct21_S10_L001_R1_001.fastq
trim_galore ZL-MeValidate-enPCR-19-Oct21_S13_L001_R1_001.fastq
trim_galore ZL-MeValidate-enPCR-19-Oct21_S11_L001_R1_001.fastq
trim_galore ZL-MeValidate-enPCR-19-Oct21_S18_L001_R1_001.fastq
trim_galore ZL-MeValidate-enPCR-19-Oct21_S22_L001_R1_001.fastq

bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/hg19/bismark ZL-MeValidate-enPCR-15-Oct21_S6_L001_R1_001.fastq_qual_trimmed.fastq

grep GCGAGAAGGCTAG *fastq


for i in `ls *fastq`
do
bismark -q --phred33-quals --bowtie2 -n 1 -l 20 /home/shg047/db/aligndb/mm9 ZL-MeValidate-enPCR-10-Oct21_S1_L001_R1_001.fastq
done



bismark 

bismark -q --phred64-quals -n 1 -l 40 --non_directional /data/genomes/homo_sapiens/GRCh37/s_1_sequence.txt


# 2015-10-26
ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/ ./

/home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/
/home/k4zhang/my_oasis_tscc/MONOD/Ecker_Tissue_WGBS/BAMfiles
/home/k4zhang/my_oasis_tscc/MONOD/N37_WGBS/BAMfiles
/home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/
/home/k4zhang/my_oasis_tscc/MONOD/tumor_WGBS/BAMfiles/
/home/k4zhang/my_oasis_tscc/MONOD/RRBS_merge/Bam_Merged_Feb2015


# 2015-10-22
# chr17:16955467-16955609 in Colon_primary_tumor.chr17.sorted.clipped.bam
samtools tview -p chr17:16955467 /home/k4zhang/my_oasis_tscc/MONOD/tumor_WGBS/BAMfiles/Colon_primary_tumor.chr17.sorted.clipped.bam  /home/shg047/db/hg19.fa
# chr11:110910907-110910951 in middle-age.chr11.rmdup.bam
samtools tview -p chr11:110910907 /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam /home/shg047/db/hg19.fa
# grep chr11:110910907-110910951 in hapInfo files of middle-age.chr11.hapInfo.txt
grep chr11:110910907-110910951 /home/kunzhang/CpgMIP/MONOD/Data/WGBS_data/mld_block_hapInfo_July2015/by_chrosomes/middle-age.chr11.hapInfo.txt
grep chr11:110910907-110910951  /home/shg047/monod/hap/wgbs/All_chromosomes_combined/WB_middle-age.all_chrs.hapInfo.txt

# RRBS
samtools tview -p chr17:47030368-47030542 /home/k4zhang/my_oasis_tscc/MONOD/RRBS_merge/Bam_Merged_Feb2015/6-P-10.merged.bam /home/shg047/db/hg19.fa
samtools tview -p chr22:37908213-37908261 /home/k4zhang/my_oasis_tscc/MONOD/RRBS_merge/Bam_Merged_Feb2015/6-P-10.merged.bam /home/shg047/db/hg19.fa
samtools tview -p chr20:44441518-44441605 /home/k4zhang/my_oasis_tscc/MONOD/RRBS_merge/Bam_Merged_Feb2015/6-P-10.merged.bam /home/shg047/db/hg19.fa

# BSPP
samtools tview -p chr1:117909025 /home/k4zhang/my_oasis_tscc/MONOD/WGBS_BSPP/BAMfiles/6P-10.sorted.clipped.rmdup.bam /home/shg047/db/hg19.fa

samtools tview -p chr1:10003741 /home/k4zhang/my_oasis_tscc/MONOD/hESC_WGBS/BAMfiles/methylC-seq_h1-npc_r1.chr1.rmdup.bam /home/shg047/db/hg19.fa


# achieve regions from methylation mhl files
cd /home/shg047/monod/dec
awk -F'[:-\t]' 'NR>1 {print $1,$2,$3}' OFS="\t" WGBS_methHap_load_matrix_20Oct2015.txt > WGBS_methHap_load_matrix_20Oct2015.bed
# transfer from genome-miner to tscc to extract aml
cd /home/shg047/monod/aml
scp shg047@genome-miner.ucsd.edu:/home/shg047/monod/dec/WGBS_methHap_load_matrix_20Oct2015.bed ./

PileOMeth extract  -q 10 -p 5 -l WGBS_methHap_load_matrix_20Oct2015.bed /home/shg047/db/hg19.fa  /home/k4zhang/my_oasis_tscc/MONOD/tumor_WGBS/BAMfiles/Colon_primary_tumor.chr17.sorted.clipped.bam  -o Colon_primary_tumor.chr17.sorted.clipped
PileOMeth extract  -q 10 -p 5 -l WGBS_methHap_load_matrix_20Oct2015.bed /home/shg047/db/hg19.fa  /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam  -o abc
PileOMeth extract  -q 10 -p 5 -r chr11:110910907-110910908 /home/shg047/db/hg19.fa /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam  -o abc
PileOMeth extract  -q 10 -p 5 -r chr11:72295556-72295573 /home/shg047/db/hg19.fa /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam  -o abc
PileOMeth extract  -q 10 -p 5 -l WGBS_methHap_load_matrix_20Oct2015.bed.head /home/shg047/db/hg19.fa  /home/k4zhang/my_oasis_tscc/MONOD/whole_blood_WGBS/BAMfiles/middle-age.chr11.rmdup.bam  -o abc
chr11:72295556-72295573

my $ref="/home/shg047/db/hg19.fa";
my $region="chr22:16050002-17050002";
my $bamfile="methylC-seq_h1_r2.chr22.rmdup.bam";
system("PileOMeth extract -r $region $ref $bam > aa.txt");


cd /home/shg047/monod/hap/wgbs/All_chromosomes_combined
screen -s a
perl /home/shg047/monod/hap/get_avemeth_matrix_16Oct2015.pl  > WGBS_aveBedMeth_load_matrix_20Oct2015

ftp address:  137.189.133.62
ftp 137.189.133.62
Username: plamethy
Password: de$*d@s3
prompt # to turn prompt off
mget -c -r *  ./


mv middle-age.all_chrs.hapInfo.txt WB_middle-age.all_chrs.hapInfo.txt
mv new-born.all_chrs.hapInfo.txt WB_new-born.all_chrs.hapInfo.txt
mv centenarian.all_chrs.hapInfo.txt WB_centenarian.all_chrs.hapInfo.txt

perl ../../get_avemeth_matrix_16Oct2015.pl > WGBS_aveMeth_load_matrix_20Oct2015.txt

wget -m --ftp-user=plamethy --ftp-password='de$*d@s3' ftp://137.189.133.62/




