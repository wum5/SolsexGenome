#!/bin/bash

#PBS -N MSK
#PBS -l nodes=1:ppn=4,walltime=12:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


## modules to load
module load java
module load trimmomatic
module load bwa
module load samtools
module load python
module load biopython

## environment setup
cd /N/dc2/projects/solanumgenome/sapp_/kmer-analyses  # the working directory
PATH=$PATH:/N/dc2/projects/solanumgenome/scripts   # the directory saving the custom python scripts
PATH=$PATH:~/softwares/seqtk   # the PATH to seqtk program
OD=/N/dc2/projects/solanumgenome/sapp_/reads/pop-dna  # the directory saving the read FASTQ files
ThreadN=4   # number of threads for parallel computing
ProcStep=1234   # step control

MSK=/N/project/solsex/_sapp/corrected_pacbio/sex_specific_reads/MSK-dna.final.wtRev.txt
FSK=/N/project/solsex/_sapp/corrected_pacbio/sex_specific_reads/FSK-dna.final.wtRev.txt

## extract PE reads with FSK
if [[ $ProcStep = *"1"* ]];
then
    mkdir FSK_pattern
    cd FSK_pattern

    for file in $OD/f_*_paired_*.fastq.gz;
    do
        zgrep -B1 -F -f ${FSK} ${file} > FSR.out

        index=$(basename $file ".fastq.gz")
        cat FSR.out | grep "^@" | cut -d" " -f1 | uniq > ../${index}.FSK.reads.txt
        rm FSR.out
    done

    cd ..

    for file in f_*_paired_R1.FSK.reads.txt;
    do
        index=$(basename $file "_paired_R1.FSK.reads.txt")
        echo ${index}
        cat ${index}_paired_R1.FSK.reads.txt ${index}_paired_R2.FSK.reads.txt | \
        sort -u > ${index}.FSK_reads.txt
        sed -i 's/@//g' ${index}.FSK_reads.txt
    done

    rm *_paired_*.FSK.reads.txt
    rm FSK_pattern/*

    mv *.FSK_reads.txt FSK_pattern/
    mkdir FSK_reads

    for file in FSK_pattern/*.FSK_reads.txt;
    do
        index=$(basename ${file} ".FSK_reads.txt")
        seqtk subseq ${OD}/${index}_paired_R1.fastq.gz ${file} > FSK_reads/${index}_R1.FSK.fastq
        seqtk subseq ${OD}/${index}_paired_R2.fastq.gz ${file} > FSK_reads/${index}_R2.FSK.fastq
    done

    cat FSK_reads/*_R1.FSK.fastq > FSK_reads.R1.fastq
    cat FSK_reads/*_R2.FSK.fastq > FSK_reads.R2.fastq
fi


## map the female specific reads back to the assembled genome sapp.v1.0
if [[ $ProcStep = *"2"* ]];
then
    mkdir fsk2assembly && cd fsk2assembly
    REF=/N/dc2/projects/solanumgenome/KEYFILE/sapp.scaffolds.fasta

    bwa index ${REF}
    bwa mem ${REF} ../FSK_reads.R1.fastq ../FSK_reads.R2.fastq \
    -t ${ThreadN} | samtools sort -@ ${ThreadN} -O BAM -o FSK_map2assembly.sort.bam -

    samtools flagstat FSK_map2assembly.sort.bam > FSK_map2assembly.summary
    samtools depth -q 20 -Q 20 -aa FSK_map2assembly.sort.bam > FSK_map2assembly.DP.txt

    window_depth.py -i FSK_map2assembly.DP.txt -w 10000 -m 300 > FSK_map2assembly.10kDP.txt
    window_depth.py -i FSK_map2assembly.DP.txt -w 1000 -m 300 > FSK_map2assembly.1kDP.txt
    cd ..
fi


## extract PE reads with MSK
if [[ $ProcStep = *"3"* ]];
then
    mkdir MSK_pattern
    cd MSK_pattern

    for file in $OD/m_*_paired_*.fastq.gz;
    do
        zgrep -B1 -F -f ${MSK} ${file} > MSR.out

        index=$(basename $file ".fastq.gz")
        cat MSR.out | grep "^@" | cut -d" " -f1 | uniq > ../${index}.MSK.reads.txt
        rm MSR.out
    done

    cd ..

    for file in m_*_paired_R1.MSK.reads.txt;
    do
        index=$(basename $file "_paired_R1.MSK.reads.txt")
        echo ${index}
        cat ${index}_paired_R1.MSK.reads.txt ${index}_paired_R2.MSK.reads.txt | \
        sort -u > ${index}.MSK_reads.txt
        sed -i 's/@//g' ${index}.MSK_reads.txt
    done

    rm *_paired_*.MSK.reads.txt
    rm MSK_pattern/*
    mv *.MSK_reads.txt MSK_pattern/

    mkdir MSK_reads

    for file in MSK_pattern/*.MSK_reads.txt;
    do
        index=$(basename ${file} ".MSK_reads.txt")
        seqtk subseq ${OD}/${index}_paired_R1.fastq.gz ${file} > MSK_reads/${index}_R1.MSK.fastq
        seqtk subseq ${OD}/${index}_paired_R2.fastq.gz ${file} > MSK_reads/${index}_R2.MSK.fastq
    done

    cat MSK_reads/*_R1.MSK.fastq > MSK_reads.R1.fastq
    cat MSK_reads/*_R2.MSK.fastq > MSK_reads.R2.fastq
fi


## map the male specific reads back to the assembled genome sapp.v1.0
if [[ $ProcStep = *"4"* ]];
then
    mkdir msk2assembly && cd msk2assembly
    REF=/N/dc2/projects/solanumgenome/KEYFILE/sapp.scaffolds.fasta

    bwa index ${REF}
    bwa mem ${REF} ../MSK_reads.R1.fastq ../MSK_reads.R2.fastq \
    -t ${ThreadN} | samtools sort -@ ${ThreadN} -O BAM -o MSK_map2assembly.sort.bam -

    samtools flagstat MSK_map2assembly.sort.bam > MSK_map2assembly.summary
    samtools depth -q 20 -Q 20 -aa MSK_map2assembly.sort.bam > MSK_map2assembly.DP.txt

    window_depth.py -i MSK_map2assembly.DP.txt -w 10000 -m 200 > MSK_map2assembly.10kDP.txt
    window_depth.py -i MSK_map2assembly.DP.txt -w 1000 -m 200 > MSK_map2assembly.1kDP.txt
    cd ..
fi