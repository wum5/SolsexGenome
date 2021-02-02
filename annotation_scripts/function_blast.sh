#!/bin/bash
set -e
set -u
set -o pipefail


function usage(){
printf "\nDescription: this script is to prepare inputs for AHRD annotation. The required softwares including \
BLAST must be added into PATH.

$(basename "$0") [-h] [-OPTIONS] -d <working dir> -s <protein fasta file> -B <blast database dir> 

where:
   -h  show this help text
   -d  the working directory for gene annotation
   -s  the input gene protein-coding sequences 
   -B  the directory saving all the blast datasets in fasta formats (e.g. uniprot_sprot.fasta)
   -n  number of cpus for parallel processing (default=24)
   -S  the steps to processed (default=123): 1: prepare; 2: Blast; 3: Interproscan\n\n"
}


######### Default parameters #########
WORKDIR=''
ProtSeq=''
BLASTDB=''
ThreadN=24
ProcStep=12


######### Parse input #########
while getopts 'hd:s:B:n:S:' option; do
  case "$option" in
    h) usage
       exit
       ;;
    d) WORKDIR=${OPTARG}
       ;;
    s) ProtSeq=${OPTARG}
       ;;
    B) BLASTDB=${OPTARG}
       ;;
    n) ThreadN=${OPTARG}
       ;;
    S) ProcStep=${OPTARG}
       ;;
    *) usage
       exit
       ;;
  esac
done
shift $((OPTIND-1))



if [ -z "${WORKDIR}" ] || [ -z "${ProtSeq}" ] || [ -z "${BLASTDB}" ] ; 
then
    usage
    exit
fi


cd $WORKDIR 



############ Prepare input sequences and databases ############
if [[ $ProcStep = *"1"* ]]; 
then
    Prefix=$(basename ${ProtSeq} .fasta)
    cp ${ProtSeq} ./${Prefix}
    splitFasta.pl -i ${Prefix} -n ${ThreadN}
    rm ${Prefix}
fi



############ BLAST homology search ############
if [[ $ProcStep = *"2"* ]]; 
then
    for DB in ${BLASTDB}/*;
    do
        DB_Prf=$(basename ${DB} .fasta)
        makeblastdb -in ${DB} -dbtype prot -out ${DB_Prf}

        find *_c*.fasta | xargs -n 1 -I{} basename {} .fasta | xargs -P ${ThreadN} -I{} \
        blastp -outfmt 6 -query {}.fasta -db ${DB_Prf} -out {}_${DB_Prf}.blastp

        cat *_${DB_Prf}.blastp > ${Prefix}_${DB_Prf}.output.blastp
    done
fi

