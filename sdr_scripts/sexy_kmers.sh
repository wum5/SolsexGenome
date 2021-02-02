#!/bin/bash

#PBS -N JELLY-DNA
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=12gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


set +o posix
module load jellyfish
module load bwa
module load samtools
module load gcc

cd /N/project/solsex/_sapp/dna-kmers

ThreadN=4
ProcStep=1234   # step control

TD=/N/project/solsex/_sapp/reads/illumina  # the two libraries (f_670-34 and m_670-1190) with high coverage
OD=/N/project/solsex/_sapp/reads/pop-data  # the rest libraries (5 males and 5 females) with low coverage


if [[ $ProcStep = *"1"* ]];
then
    jellyfish count <(zcat $TD/f_670-34_R1.fastq.gz) <(zcat $TD/f_670-34_R2.fastq.gz) \
    -m 30 --min-qual-char=? -t ${ThreadN} -C -s 2G -o f_670-34.jf
    jellyfish histo -t ${ThreadN} f_670-34.jf > f_670-34.histo

    jellyfish count <(zcat $TD/m_670-1190_R1.fastq.gz) <(zcat $TD/m_670-1190_R2.fastq.gz) \
    -m 30 --min-qual-char=? -t ${ThreadN} -C -s 2G -o m_670-1190.jf
    jellyfish histo -t ${ThreadN} m_670-1190.jf > m_670-1190.histo

    jellyfish dump m_670-1190.jf -c | sort -k1,1 > m_670-1190.kmers.txt
    jellyfish dump f_670-34.jf -c | sort -k1,1 > f_670-34.kmers.txt
    join -a1 -a2 -o 0 1.2 2.2 -e "0" f_670-34.kmers.txt m_670-1190.kmers.txt > merged.kmers.txt
fi


if [[ $ProcStep = *"2"* ]];
then
    for file in ${OD}/*_R1.fq;
    do
        #ID=$(basename ${file} "_R1.fastq.gz")
        #jellyfish count <(zcat ${PD}/${ID}_R1.fastq.gz) <(zcat ${PD}/${ID}_R2.fastq.gz) \
        #-m 30 --min-qual-char=? -t ${ThreadN} -C -s 2G -o ${ID}.jf

        ID=$(basename ${file} "_R1.fq")
        jellyfish count ${OD}/${ID}_R1.fq ${OD}/${ID}_R2.fq \
        -m 30 --min-qual-char=? -t ${ThreadN} -C -s 2G -o ${ID}.jf

        #ID=$(basename ${file} ".jf")
        #jellyfish dump ${PD}/${ID}.jf -c | sort -k1,1 > ${ID}.kmers.txt

        #join -a1 -a2 -o 0 1.2 2.2 -e "0" malpool.kmers.txt fempool.kmers.txt | \
        #join -a1 -a2 -o 0 1.2 2.2 2.3 -e "0" herpool.kmers.txt - > merged.kmers.txt
    done
fi


if [[ $ProcStep = *"3"* ]];
then
    awk '$2>=5 && $3==0' merged.kmers.txt > FSK.raw.txt
    awk '$2==0 && $3>=5' merged.kmers.txt > MSK.raw.txt

    join -a1 -a2 -o 0 1.2 2.2 -e "0" f_716-25.kmers.txt f_716-4.kmers.txt | \
    join -a1 -a2 -o 0 1.2 2.2 2.3 -e "0" f_716-23.kmers.txt - | \
    join -a1 -a2 -o 0 1.2 2.2 2.3 2.4 -e "0" f_670-72.kmers.txt - | \
    join -a1 -a2 -o 0 1.2 2.2 2.3 2.4 2.5 -e "0" f_670-53.kmers.txt - > f_pop.kmers.txt

    join -a1 -a2 -o 0 1.2 2.2 -e "0" m_716-27.kmers.txt m_716-6.kmers.txt | \
    join -a1 -a2 -o 0 1.2 2.2 2.3 -e "0" m_716-24.kmers.txt - | \
    join -a1 -a2 -o 0 1.2 2.2 2.3 2.4 -e "0" m_670-59.kmers.txt - | \
    join -a1 -a2 -o 0 1.2 2.2 2.3 2.4 2.5 -e "0" m_670-105.kmers.txt - > m_pop.kmers.txt

    join -a1 -a2 -o 0 1.2 1.3 1.4 1.5 1.6 2.2 2.3 2.4 2.5 2.6 -e "0" \
    f_pop.kmers.txt m_pop.kmers.txt > all_pop.kmers.txt

    awk '$7==0 && $8==0 && $9==0 && $10==0 && $11==0' all_pop.kmers.txt > f_pop.filter.kmers
    awk '$2==0 && $3==0 && $4==0 && $5==0 && $6==0' all_pop.kmers.txt > m_pop.filter.kmers

    join -o 1.1 1.2 2.2 2.3 2.4 2.5 2.6 FSK.raw.txt f_pop.filter.kmers > FSK.filtered.txt
    join -o 1.1 1.3 2.7 2.8 2.9 2.10 2.11 MSK.raw.txt m_pop.filter.kmers > MSK.filtered.txt
fi


if [[ $ProcStep = *"4"* ]];
then
    awk 'BEGIN{COUNT=1} $3>0 && $4>0 && $5>0 && $6>0 && $7>0 \
    {print ">seq"COUNT"\n"$1; COUNT=COUNT+1}' FSK.filtered.txt > FSK-dna.final.fasta

    bioawk -c fastx '{print ">"$name; print $seq; print">"$name"_rv"; print revcomp($seq)}' \
    FSK-dna.final.fasta > FSK-dna.final.wtRev.fasta

    grep -v '^>' FSK-dna.final.wtRev.fasta > FSK-dna.final.wtRev.txt

    awk 'BEGIN{COUNT=1} $3>0 && $4>0 && $5>0 && $6>0 && $7>0 \
    {print ">seq"COUNT"\n"$1; COUNT=COUNT+1}' MSK.filtered.txt > MSK-dna.final.fasta

    bioawk -c fastx '{print ">"$name; print $seq; print">"$name"_rv"; print revcomp($seq)}' \
    MSK-dna.final.fasta > MSK-dna.final.wtRev.fasta

    grep -v '^>' MSK-dna.final.wtRev.fasta > MSK-dna.final.wtRev.txt
fi