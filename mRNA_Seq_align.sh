#!/bin/bash
export $PATH

# for i in *.bt2; do echo $i; done

for i in *.fastq;

do

#mkdir ./${i%.fastq}

tophat2 -G Bdistachyon_314_v3.1.gene.gff3 -p 4 Bdistachyon_314_v3.0 $i ;

# --transcriptome-only

# cufflinks [options] <hits.sam> or .bam

cufflinks -o ${i%.fastq}_cl ${i%.fastq}/accepted_hits.bam

done

# output into $i_cl is going to be a .gtf file , take all .gtfs and put into assemblies file.

#cuffmerge [Options] <assembly_GTF_list.txt>

#cuffmerge ./assemblies.txt

# cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam>

#cuffdiff -L WT,mutant1  ./ASM/outputfromcuffmerge.gtf ~/Bd_5GS_RNAseq/0_2_R1/accepted_hits.bam,~/Bd_5GS_RNAseq/0_2_R2/accepted_hits.bam mutrep1.bam,mutrep2.bam

# gene_exp.diff   isoforms.diff
