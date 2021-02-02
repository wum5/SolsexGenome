g#!/bin/bash

#PBS -N Variant
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=32gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


### set up environment
module load picard/2.14.0 
module load samtools/1.5
module load vcftools
PICARD_PATH=/N/soft/rhel7/picard/2.14.0
GATK_PATH=/N/u/wum5/Carbonate/softwares/GenomeAnalysisTK-3.8
ANGSD_PATH=/N/u/wum5/Carbonate/softwares/angsd
GENGER_PATH=/N/u/wum5/Carbonate/softwares/genomics_general  # https://github.com/simonhmartin/genomics_general


### set up input files or parameters
GENOME=/N/dc2/projects/solanumgenome/KEYFILE/sapp.scaffolds.fasta
JAVA_OPT="-Xmx8G"
ThreadN=4
ProcStep=123456   # step control


### move to working directory
cd /N/dc2/projects/solanumgenome/sapp/dna-bams   # the direcotry having the BAM files of DNA alignment of all plant individuals


### remove the duplicate sequences
if [[ $ProcStep = *"1"* ]];
then
    for file in *.sort.bam;
    do
        index=$(basename $file .sort.bam)

        java ${JAVA_OPT} -jar ${PICARD_PATH}/picard.jar MarkDuplicates \
        TMP_DIR=tmp I=${index}.sort.bam O=${index}.dedup.bam \
        METRICS_FILE=${index}.dedup.metrics.txt ASSUME_SORTED=true \
        REMOVE_DUPLICATES=false TAGGING_POLICY=All

        mv *.dedup.metrics.txt tmp/
    done
fi


### index the genome
if [[ $ProcStep = *"2"* ]];
then
    cp ${GENOME} REF.fasta
    samtools faidx REF.fasta

    java ${JAVA_OPT} -jar ${PICARD_PATH}/picard.jar CreateSequenceDictionary \
    R=REF.fasta O=REF.dict


    ### add header to bam files
    index=0
    for file in *.dedup.bam;
    do
        index=$((index+1))
        sample=$(basename $file .dedup.bam)
        java ${JAVA_OPT} -jar ${PICARD_PATH}/picard.jar AddOrReplaceReadGroups \
        I=${sample}.dedup.bam O=${sample}.dedup.renamed.bam \
        RGID=${index} RGLB="lib"${index} RGPL=illumina RGPU="unit"${index} RGSM=${sample}
    done


    ### index the renamed bam files
    for file in *.dedup.renamed.bam;
    do
        java ${JAVA_OPT} -jar ${PICARD_PATH}/picard.jar BuildBamIndex I=${file}
    done
fi

    
### re-aligner target creator
if [[ $ProcStep = *"3"* ]];
then
    INPUT=$(find *.dedup.renamed.bam | xargs -n 1 -I{} echo "-I {}")

    java ${JAVA_OPT} -jar ${GATK_PATH}/GenomeAnalysisTK.jar -T RealignerTargetCreator \
    -R REF.fasta ${INPUT} -o all_samples.intervals -nct ${ThreadN}


    ### indel realignment
    for file in *.dedup.renamed.bam;
    do
        index=$(basename $file .dedup.renamed.bam)

        java ${JAVA_OPT} -jar ${GATK_PATH}/GenomeAnalysisTK.jar -T IndelRealigner \
        -R REF.fasta -I ${index}.dedup.renamed.bam \
        -targetIntervals all_samples.intervals -o ${index}.realigned.bam
    done
fi


### calling variants using ANGSD with GATK genotype likelihood model
if [[ $ProcStep = *"4"* ]];
then
    ${ANGSD_PATH}/angsd -GL 2 -nThreads ${ThreadN} -doGlf 2 \
    -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam bam.filelist \
    -minMapQ 30 -minQ 20 -minMaf 0.05 -minInd 8 -out genolike
fi


### calling variants for each sample
if [[ $ProcStep = *"5"* ]];
then
    for file in *.realigned.bam;
    do
        index=$(basename $file .realigned.bam)
    
        java ${JAVA_OPT} -jar ${GATK_PATH}/GenomeAnalysisTK.jar -T HaplotypeCaller \
        -R REF.fasta -I ${index}.realigned.bam -o ${index}.raw.g.vcf \
        --emitRefConfidence GVCF -nct ${ThreadN}
    done
fi


### joint genotyping across samples
if [[ $ProcStep = *"6"* ]];
then
    INPUT=$(find *.raw.g.vcf | xargs -n 1 -I{} echo "-V {}")

    java ${JAVA_OPT} -jar ${GATK_PATH}/GenomeAnalysisTK.jar -T GenotypeGVCFs \
    -allSites -R REF.fasta ${INPUT} -o allsamples.allSites.vcf

    python ${GENGER_PATH}/VCF_processing/parseVCF.py -i allsamples.allSites.vcf --skipIndels \
    --minQual 30 --gtf flag=DP min=5 | gzip > allsamples.HD5.geno.gz

    zgrep -v 'N/N' allsamples.HD5.geno.gz | gzip > allsamples.HD5-filtered.geno.gz

    python ${GENGER_PATH}/popgenWindows.py -w 10000 -m 1000 \
    -g allsamples.HD5-filtered.geno.gz \
    -o 670_vs_716.HD5.csv.gz -f phased -T 4 \
    -p popA f_670-34,f_670-53,f_670-72,m_670-105,m_670-1190,m_670-59 \
    -p popB f_716-23,f_716-25,f_716-4,m_716-24,m_716-27,m_716-6

    python ${GENGER_PATH}/popgenWindows.py -w 10000 -m 1000 \
    -g allsamples.HD5-filtered.geno.gz \
    -o m_pop_diff.HD5.csv.gz -f phased -T 4 \
    -p popA m_670-105,m_670-1190,m_670-59 \
    -p popB m_716-24,m_716-27,m_716-6

    python ${GENGER_PATH}/popgenWindows.py -w 10000 -m 1000 \
    -g allsamples.HD5-filtered.geno.gz \
    -o f_pop_diff.HD5.csv.gz -f phased -T 4 \
    -p popA f_670-34,f_670-53,f_670-72 \
    -p popB f_716-23,f_716-25,f_716-4

    python ${GENGER_PATH}/popgenWindows.py -w 10000 -m 1000 \
    -g allsamples.HD5-filtered.geno.gz \
    -o 670_sex-diff.HD5.csv.gz -f phased -T 4 \
    -p popA f_670-34,f_670-53,f_670-72 \
    -p popB m_670-105,m_670-1190,m_670-59

    python ${GENGER_PATH}/popgenWindows.py -w 10000 -m 1000 \
    -g allsamples.HD5-filtered.geno.gz \
    -o 716_sex-diff.HD5.csv.gz -f phased -T 4 \
    -p popA f_716-23,f_716-25,f_716-4 \
    -p popB m_716-24,m_716-27,m_716-6
fi
