#! bin/bash
# Script is from Dr. Weilin Pu and the input configure file is as the following:

#! bin/bash

declare -a GroupName
declare -a GroupSamplesName

i=1
IFS=""

for LINE in `cat ./sampleInfo.txt`
do
	split1=`echo "${LINE}" | cut -d " " -f1` 
	split2=`echo "${LINE}" | cut -d " " -f2`
	GroupName[$i]="$split1"
	GroupSamplesName[$i]="$split2"
	((i++))
done

Ngroups=${#GroupName[*]}



for ((j=1; j<=${#GroupName[*]}; j++))
do
	echo ${GroupName[$j]}
	echo ${GroupSamplesName[$j]}
done


tophat2 -p 8 -G ./genes.gtf --transcriptome-index=transcriptome_data/known ./genome


for i in `ls | grep R1.fastq.gz`
do
newName="`echo $i | sed -n 's/_R1.fastq.gz//p'`"
file1="${newName}_R1.fastq.gz"
file2="${newName}_R2.fastq.gz"
output1="${newName}_thout"
output2="${newName}_clout"
output3="${output1}/accepted_hits.bam"
output4="${output2}/transcripts.gtf"
tophat2 -p 8 -o $output1 --transcriptome-index=transcriptome_data/known ./genome $file1 $file2
cufflinks -p 8 -u -b genome.fa -G genes.gtf -o $output2 $output3
echo ${output4} >>test_assemblies.txt
done


cuffmerge -g genes.gtf -s genome.fa -p 8 test_assemblies.txt

for ((c=1; c<=$Ngroups; c++))
do
for ((d=$[c+1]; d<=$Ngroups; d++))
do

#Get the samples of the first group 
i=1  
while((1==1))  
do  
        split=`echo "${GroupSamplesName[$c]}" | cut -d ";" -f$i`  
        if [ "$split" != "" ]  
        then  
                
                samples1[$i]="$split"  
		((i++))  
        else  
               break  
        fi  
done  



#Get the samples of the second group
i=1  
while((1==1))  
do  
        split=`echo "${GroupSamplesName[$d]}" | cut -d ";" -f$i`  
        if [ "$split" != "" ]  
        then  
                
                samples2[$i]="$split"  
		((i++))  
        else  
               break  
        fi  
done  


#Get the combined input of the first group
testoutput1=${samples1[1]}"_thout/accepted_hits.bam"
for ((j=2; j<= ${#samples1[*]}; j++ ))
do
testoutput1=${testoutput1}",""${samples1[$j]}_thout/accepted_hits.bam"
done


#Get the combined input of the second group 
testoutput2=${samples2[1]}"_thout/accepted_hits.bam"
for ((j=2; j<= ${#samples2[*]}; j++ ))
do
testoutput2=${testoutput2}",""${samples2[$j]}_thout/accepted_hits.bam"
done

diffname="${GroupName[$c]}_${GroupName[$d]}_diff_out"
cuffdiff -o $diffname -b genome.fa -p 8 -L "${GroupName[$c]},${GroupName[$d]}" -u merged_asm/merged.gtf $testoutput1 $testoutput2

done
done

