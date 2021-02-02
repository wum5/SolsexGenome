#!/bin/bash
set -e
set -u
set -o pipefail


function usage(){
printf "\nDescription: this script is to do quality control of base calling in the genome assembly.
Make sure the four programs (BWA, SAMTOOLS, ANGSD, and REFEREE) are installed and added to PATH.
Three are 3 steps in this pipeline: 1) map reads using BWA and SAMtools;
2) calculate genotype likelihoods using ANGSD; 3) calculate base quality score.
$(basename "$0") [-h] [-OPTIONS] -d <working dir> -g <genome file>  \
-f <foward reads> -r <reverse reads>
where:
   -h  show this help text
   -d  the working directory
   -g  genome assembly FASTA file
   -f  forward reads FASTQ file
   -r  reverse reads FASTQ file
   -t  number of threads for parallel processing (default=16)
   -S  the steps to processed (default=123)\n\n"
}


######### Default parameters #########
WORKDIR=''
GENOME=''
READ1=''
READ2=''
ThreadN=4
ProcStep=123


######### Parse input #########
while getopts 'h:d:g:f:r:t:S:' option; do
  case "$option" in
    h) usage
       exit
       ;;
    d) WORKDIR=${OPTARG}
       ;;
    g) GENOME=${OPTARG}
       ;;
    f) READ1=${OPTARG}
       ;;
    r) READ2=${OPTARG}
       ;;
    t) ThreadN=${OPTARG}
       ;;
    S) ProcStep=${OPTARG}
       ;;
    *) usage
       exit
       ;;
  esac
done
shift $((OPTIND-1))



if [ -z "${WORKDIR}" ] || [ -z "${GENOME}" ] || [ -z "${READ1}" ] || [ -z "${READ2}" ];
then
    usage
    exit
fi

cd ${WORKDIR}
INDEX=$(basename ${READ1} .fastq.gz)


### DNA mapping
if [[ ${ProcStep} = *"1"* ]];
then
    if [ ! -f REF ]; then
        cp ${GENOME} REF
        bwa index REF
        samtools faidx REF
    fi

    bwa mem REF ${READ1} ${READ2} -t ${ThreadN} | \
    samtools sort -@ ${ThreadN} -O BAM -o ${INDEX}.sort.bam -
fi


### Calculating genotype log-likelihoods with ANGSD
if [[ ${ProcStep} = *"2"* ]];
then
    angsd -GL 2 -i ${INDEX}.sort.bam -ref REF -minQ 0 -doGlf 4 -out ${INDEX}-gl
fi


### Calculate reference quality scores with Referee
if [[ ${ProcStep} = *"3"* ]];
then
    referee.py -gl ${INDEX}-gl.glf.gz -ref REF --fastq -o ${INDEX}-ref
fi


