#!/bin/bash

#PBS -N RNA-MF
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=4gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module load trimmomatic
module load hisat2
module load samtools
module load stringtie
module load htseq


### Sapp input
ID=/N/dc2/projects/solanumgenome/sapp_/reads/all-rna  # directory saving input fastq files
OD=/N/dc2/projects/solanumgenome/sapp_/rna-bams  # directory saving output bam files
SD=/N/dc2/projects/solanumgenome/scripts  # directory saving scripts
GENOME=/N/dc2/projects/solanumgenome/KEYFILE/sapp.scaffolds.fasta  # genome assembly sequence
ThreadN=8


### RNA mapping 
cd ${OD}
hisat2-build -p ${ThreadN} ${GENOME} REF

for file in ${ID}/*_R1.fastq.gz;
do
    index=$(basename $file _R1.fastq.gz)

    mkdir tmp
    cd tmp

    trimmomatic-0.36 PE -phred33 -threads ${ThreadN} \
    ${ID}/${index}_R1.fastq.gz ${ID}/${index}_R2.fastq.gz \
    ${index}_paired_R1.fastq.gz ${index}_unpaired_R1.fastq.gz \
    ${index}_paired_R2.fastq.gz ${index}_unpaired_R2.fastq.gz\
    ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

    cd ..

    hisat2 -p ${ThreadN} --dta -x REF \
    -1 tmp/${index}_paired_R1.fastq.gz \
    -2 tmp/${index}_paired_R2.fastq.gz \
    -S ${index}.sam

    samtools sort -@ ${ThreadN} -o ${index}.bam ${index}.sam

    stringtie -p ${ThreadN} -o ${index}.gtf ${index}.bam
done


## Extract transcript sequences into a FASTA file used for MAKER gene prediction
stringtie --merge -p ${ThreadN} -o stringtie_merged.gtf *

~/softwares/gffread-0.9.12/gffread -w stringtie.transcripts.fasta -g ${GENOME} stringtie_merged.gtf

