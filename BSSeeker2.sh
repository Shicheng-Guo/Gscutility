#! /bin/bash
#=============Set Working Directory====================
#=============Get the file name========================
Folder="/home/shg047/oasis/Minghua2016/fastq"
cd $Folder
Output_file="Sample.txt"
:> $Output_file
for file_a in ${Folder}/*fastq;
do
    temp_file=`basename $file_a`
    echo $temp_file>>$Output_file
done

#============Start the reference build =================
python bs_seeker2-build.py -f /home/shg047/db/hg19/hg19.fa --aligner=bowtie2
#============Start the alignment  and methylation ratio calculation process ================
declare -a lines
lines=($(cat Sample.txt))
N_sample=$((${#lines[@]}/2))

for((i=0;i<${N_sample};i++));do
     idx1=$((2*$i))
     idx2=$((2*$i+1))
     file1=${lines[$idx1]}
     file2=${lines[$idx2]}
     SampleID1=${file1%_R*}
     SampleID2=${file1%.fa*}
     outputname=${SampleID1}".bam"
     python /home/shg047/software/BSseeker2/bs_seeker2-align.py -t Y -1 ${file1} -2 ${file2} --aligner=bowtie2 -o ${outputname} -f bam -g /home/shg047/db/hg19/hg19.fa 1>&2 2> ${SampleID1}.log
	 python /home/shg047/software/BSseeker2/bs_seeker2-call_methylation.py -i ${outputname} -o ${SampleID1} --txt --db /home/shg047/software/BSseeker2/bs_utils/reference_genomes/hg19.fa_bowtie2
	 echo $outputname
done

#==========Start to do analysis based on the CGmap data generated with R
Rscript /home/puweilin/Software/MethylRatio.R

#==========Finished The Methylation Sequencing and Analysis Proecss===============================================





