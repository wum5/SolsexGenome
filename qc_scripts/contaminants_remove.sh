#!/bin/bash
set -e
set -u
set -o pipefail


function usage(){
printf "\nThis script is to remove the potential contaminant scaffolds from the genome assembly, \
by searching the databases from human, bacteria and virus but retaining the scaffolds that have \
high similarity to the closely related species. You must ensure the Deconseq program has been \
correctly installed and set in the paths file, and the corresponding databases have been downloaded \
before running the script.
$(basename "$0") [-h] [-OPTIONS] -g <genome file> -r <retention database> -R human,bacteria,viral
where:
    -h  show this help text
    -d  PATH to the working directory
    -g  the input genome in fasta file
    -r  the retention dabase
    -R  the contaminant databse (default=human,bacteria,viral)
    -n  number of cpus for parallel processing (default=16)
    -i  the alignment identify percent cutoff (default=80)
    -c  the alignment coverage percent cutoff (default=50)\n\n"
}


######### Default parameters #########
WORKDIR=''
GENOME=''
dbsRetain=''
dbsRemove=human,bacteria,viral  
ThreadN=16 
alignidt=80	
aligncov=50



######### Parse input #########
while getopts ':h:d:g:R:r:n:i:c' option; do
  case "$option" in
    h) usage
       exit
       ;;
    d) WORKDIR=${OPTARG}
       ;;
    g) GENOME=${OPTARG}
       ;;
    r) dbsRetain=${OPTARG}
       ;;
    R) dbsRemove=${OPTARG}
       ;;
    n) ThreadN=${OPTARG}
       ;;
    i) alignidt=${OPTARG}
       ;;
    c) aligncov=${OPTARG}
       ;;
    *) usage
       exit
       ;;
  esac
done
shift $((OPTIND-1))


if [ -z "${WORKDIR}" ] || [ -z "${GENOME}" ] || [ -z "${dbsRetain}" ] ; then
    usage
    exit
fi


if [ -z "${DeconSeq+x}" ] ; then
	echo -e '\nSome required program(s) have not yet been set in paths\n'
	exit
fi


exit

########## run DeconSeq to remove contaminants ########## 
cd $WORKDIR
if [ ! -d contaminants ]; then mkdir contaminants; fi
cd contaminants

genomeID=$(basename $GENOME ".fa") 
cp $GENOME $genomeID
splitFasta.pl -i $genomeID -n 32

function run_deconseq(){
dir=$(basename ${1} .fasta)  # {1}=input file
mkdir $dir
mv ${1} $dir
cd $dir  
${DeconSeq}/deconseq.pl -f ${1} -dbs ${dbsRemove} -dbs_retain ${dbsRetain} \
-i ${alignidt} -c ${aligncov} -out_dir .
cd ..
}

export -f run_deconseq
find *fasta | xargs -n 1 | xargs -P $ThreadN -I{} bash -c 'run_deconseq {}' 

cat ${genomeID}_c*/*_clean.fa > ${genomeID}"_clean.fa"
cat ${genomeID}_c*/*_cont.fa > ${genomeID}"_cont.fa"
cat ${genomeID}_c*/*_both.fa > ${genomeID}"_both.fa" 


