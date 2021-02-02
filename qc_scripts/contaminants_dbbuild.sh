#!/bin/bash
set -e
set -u
set -o pipefail


function usage(){
printf "\nThis script is to set the databases for contaminant scaffold searching. \
The contaminant datasets (bacteria, virus and human genome) are downloaded from NCBI \
and constructed into database. You must ensure 'Prinseq' and 'Deconseq' program has been \
correctly installed and set in the paths file. \
Here we used Solanum genomes as retention database. You need to \
add the module 'retain_db' in the code to download your specific retention database. \
\
\n\n\e[31mAfter running this script, you need to edit the DeconSeqConfig.pm in the 'deconseq' \
software folder (example below):\e[39m \
\nuse constant DB_DIR => '/N/dc2/projects/solanumgenome/library/deconseq_DB/'; \
\nuse constant TMP_DIR => '/N/dc2/projects/solanumgenome/sapp/drafts/contaminants_removed/tmp/'; \
\nuse constant OUTPUT_DIR => '/N/dc2/projects/solanumgenome/sapp/drafts/contaminants_removed/output/'; \
\nuse constant PROG_NAME => 'bwa64';  # should be either bwa64 or bwaMAC (based on your system architecture) \
\nuse constant PROG_DIR => '';   # should be the location of the PROG_NAME file ('' if you add it to PATH) \
\nuse constant DBS => {human => {name => 'Human Reference GRCh38', \
\n                               db => 'human'}, \
\n                     bact => {name => 'Bacterial genomes', \
\n                              db => 'bacteria_c1,bacteria_c2,bacteria_c3,bacteria_c4,bacteria_c5,bacteria_c6,bacteria_c7,bacteria_c8,bacteria_c9,bacteria_c10'},\n \
\n                     vir => {name => 'Viral genomes', \
\n                             db => 'viral'}, \
\n                     solanum => {name => 'Solanum genomes', \
\n                             db => 'lyco,tube,penn'}}; \
\nuse constant DB_DEFAULT => 'human';
$(basename "$0") [-h] [-OPTIONS] -d <working dir> 
where:
    -h  show this help text
    -d  PATH to the directory for preparing searching datasets
    -n  number of CPUs for parallel computing (default=16)\n\n"
}


######### Default parameters #########
LIBDIR=''
ThreadN=16  


######### Parse input #########
while getopts ':h:d:n' option; do
  case "$option" in
    h) usage
       exit
       ;;
    d) LIBDIR=${OPTARG}
       ;;
    n) ThreadN=${OPTARG}
       ;;
    *) usage
       exit
       ;;
  esac
done
shift $((OPTIND-1))


if [ -z "${LIBDIR}" ] ; then
    usage
    exit
fi


if [ -z "${DeconSeq+x}" ] || [ -z "${PrinSeq+x}" ] ; then
	echo -e '\nSome required program(s) have not yet been set in paths\n'
	exit
fi



########## Contaminant database ##########
# human genome
function contam_db(){
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
mv GRCh38_latest_genomic.fna human.fasta
sample_names+="human "

## a function to download genomic files
function download_genome(){
filename=$(basename ${1})'_genomic.fna.gz'
wget -P ${2} ${1}'/'$filename
}
export -f download_genome

## viral genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt
awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > viral_genomes.txt
rm assembly_summary.txt
mkdir GbVir
cat viral_genomes.txt | xargs -n 1 | xargs -P $ThreadN -I{} bash -c 'download_genome {} GbVir'
find GbVir/*.gz | xargs -n 1 | xargs -P $ThreadN -I{} gunzip {}
cat GbVir/*.fna > viral.fasta
sample_names+="viral "

## bacterial genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > bacteria_genomes.txt
rm assembly_summary.txt
mkdir GbBac
cat bacteria_genomes.txt | xargs -n 1 | xargs -P $ThreadN -I{} bash -c 'download_genome {} GbBac'
find GbBac/*.gz | xargs -n 1 | xargs -P $ThreadN -I{} gunzip {}
cat GbBac/*.fna > bacteria.fasta
mv bacteria.fasta bacteria
splitFasta.pl -verbose -i bacteria -n 10
sample_names+="bacteria "
}


########## Retention database (this is for Solanum species; change this to your system) ##########
function retain_db(){
## S.lycopersicum
wget ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/annotation/ITAG2.4_release/S_lycopersicum_chromosomes.2.50.fa.gz
gunzip S_lycopersicum_chromosomes.2.50.fa.gz 
mv S_lycopersicum_chromosomes.2.50.fa lyco.fasta
sample_names+="lyco "

## S.pennellii
wget ftp://ftp.solgenomics.net/genomes/Solanum_pennellii/penn.fasta
sample_names+="penn "

## S.tuberosum
wget ftp://ftp.solgenomics.net/genomes/Solanum_tuberosum/assembly/PGSC_DM_v4.03/PGSC_DM_v4.03_pseudomolecules.fasta.zip
wget ftp://ftp.solgenomics.net/genomes/Solanum_tuberosum/assembly/PGSC_DM_v4.03/PGSC_DM_v4.03_unanchored_regions_chr00.fasta.zip
unzip PGSC_DM_v4.03_pseudomolecules.fasta.zip
unzip PGSC_DM_v4.03_unanchored_regions_chr00.fasta.zip
cat PGSC_DM_v4.03_unanchored_regions_chr00.fasta PGSC_DM_v4.03_pseudomolecules.fasta > tube.fasta
sample_names+="tube "
}


########## Prepare the searching databases ########## 
cd $LIBDIR
if [ ! -d deconseq_DB ]; then mkdir deconseq_DB; fi
cd deconseq_DB
sample_names=()

contam_db  # create contaminant databases
retain_db  # create retention databases
echo $sample_names

function seq_split(){
cat ${1}'.fasta' | perl -p -e 's/N\n/N/' | \
perl -p -e 's/^N+//;s/N+$//;s/N{200,}/\n>split\n/' > ${1}'_split.fa'
}
export -f seq_split

if [ "$sample_names" != '' ]
then
echo $sample_names | xargs -n 1 | xargs -P $ThreadN -I{} bash -c 'seq_split {}'

echo $sample_names | xargs -n 1 | xargs -P $ThreadN -I{} -I{} ${PrinSeq}/prinseq-lite.pl \
-log -verbose -fasta {}'_split.fa' -min_len 200 -ns_max_p 10 -derep 12345 \
-out_good {}'_split_prinseq' -seq_id {} -rm_header -out_bad null

echo $sample_names | xargs -n 1 | xargs -P $ThreadN -I{} ${DeconSeq}/bwa64 index -p {} \
-a bwtsw {}'_split_prinseq.fasta'



