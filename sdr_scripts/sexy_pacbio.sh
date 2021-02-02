#!/bin/bash

#PBS -N PAC-SSR
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu



cd /N/project/solsex/_sapp/corrected_pacbio/sex_specfic_reads

# get the corrected pacbio reads from MaSurCA outputs
cp ../assembly_male/mr.41.15.17.0.029.1.fa corrected_pacbio.m_670-1190.fa
cp ../assembly_female/mr.41.15.17.0.029.1.fa corrected_pacbio.f_670-34.fa

# get the sex-specific reads based on MSK/FSK
grep -B1 -F -f MSK-dna.final.wtRev.txt corrected_pacbio.m_670-1190.fa > pacbio_MSR.fasta
grep -B1 -F -f FSK-dna.final.wtRev.txt corrected_pacbio.f_670-34.fa > pacbio_FSR.fasta

# map those sex-specific reads to the genome assembly
~/softwares/minimap2/minimap2 -x map-pb --secondary=no /N/project/solsex/_keyfiles/sapp.scaffolds.fasta pacbio_MSR.fasta > pacbio_MSR_aln.sam
~/softwares/minimap2/minimap2 -x map-pb --secondary=no /N/project/solsex/_keyfiles/sapp.scaffolds.fasta pacbio_FSR.fasta > pacbio_FSR_aln.sam

# count the read number on each mapped scaffolds
cut -f6 pacbio_MSR_aln.sam | grep 'scf' | sort | uniq -c | sort -nr > pacbio_MSR_scaffolds.txt
cut -f6 pacbio_FSR_aln.sam | grep 'scf' | sort | uniq -c | sort -nr > pacbio_FSR_scaffolds.txt

# assemble sex-specific pacbio reads
~/softwares/smartdenovo/smartdenovo.pl -p MSR -c 1 pacbio_MSR.fasta > MSR.mak
~/softwares/smartdenovo/smartdenovo.pl -p FSR -c 1 pacbio_FSR.fasta > FSR.mak

make -f MSR.mak
make -f FSR.mak

# map those sex-specific contigs to the genome assembly
~/softwares/minimap2/minimap2 -x map-pb --secondary=no /N/project/solsex/_keyfiles/sapp.scaffolds.fasta MSR.dmo.cns > MSR.dmo.cns.sam
~/softwares/minimap2/minimap2 -x map-pb --secondary=no /N/project/solsex/_keyfiles/sapp.scaffolds.fasta FSR.dmo.cns > FSR.dmo.cns.sam

