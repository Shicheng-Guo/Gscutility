# from Weilin

#! /bin/bash

#=============Set Working Directory====================
cd /home/puweilin/Software/BSseeker2-master/

#=============Get the file name========================
Folder="/home/puweilin/Software/BSseeker2-master/RawData"
Output_file="Sample.txt"
:> $Output_file
for file_a in ${Folder}/*;do
    temp_file=`basename $file_a`
    echo $temp_file>>$Output_file
done

#============Start the reference build =================

/home/puweilin/anaconda2/bin/python bs_seeker2-build.py -f /home/puweilin/Software/bowtie2-2.2.9/ChrM/Homo_sapiens.GRCh38.dna.chromosome.MT.fa --aligner=bowtie2


#============Start the alignment  and methylation ratio calculation process ================
declare -a lines
lines=($(cat Sample.txt))
N_sample=$((${#lines[@]}/2))

for((i=0;i<${N_sample};i++));do
     idx1=$((2*$i))
     idx2=$((2*$i+1))
     file1="RawData/"${lines[$idx1]}
     file2="RawData/"${lines[$idx2]}
     SampleID1=${file1%_R*}
     SampleID2=${file1%.fa*}
     outputname=${SampleID1}".bam"
    /home/puweilin/anaconda2/bin/python bs_seeker2-align.py -1 ${file1} -2 ${file2} --aligner=bowtie2 -o ${outputname} -f bam -g /home/puweilin/Software/bowtie2-2.2.9/ChrM/Homo_sapiens.GRCh38.dna.chromosome.MT.fa

    /home/puweilin/anaconda2/bin/python bs_seeker2-call_methylation.py -i ${outputname} -o ${SampleID1} --txt --db  /home/puweilin/Software/BSseeker2-master/bs_utils/reference_genomes/Homo_sapiens.GRCh38.dna.chromosome.MT.fa_bowtie2
done


#==========Start to do analysis based on the CGmap data generated with R
Rscript /home/puweilin/Software/MethylRatio.R

#==========Finished The Methylation Sequencing and Analysis Proecss===============================================





