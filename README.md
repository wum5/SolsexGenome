Solanum appendiculatum Genome Project
=================


### Overview
The scripts in this repository are to generate the outputs in Solanum appendiculatum genome project


### Contributors 
* Meng Wu
* Rafael Guerrero


### Table of Contents
* [De novo Asseemble Genome](#de-novo-assemble-genome)
* [Estimate Genome Size](#estimate-genome-size)
* [Examine Genome Quality](#examine-genome-quality)
* [Genome Annotation](#genome-annotation)
* [Gene Family Analysis](#gene-family-analysis)
* [Transcriptome Analysis](#transcriptome-analysis)
* [Sex Determination Analysis](#sex-determination-analysis)


### De novo Assemble Genome

#### Generate genome assembly using both female and male reads (the assembly descibed in manuscript)  
```
sh assembly_scripts/pool_assembly.sh 
```
Notes: You need to change the PATH in the corresponding config file


### Estimate Genome Size

#### Generate a kmer histo file using Jellyfish
```
jellyfish count -m 25 -t 8 -C -s 4G -o male_reads.counts \
<(zcat <male_reads1.fq>) <(zcat <male_reads2.fq>)

jellyfish count -m 25 -t 8 -C -s 4G -o fema_reads.counts \
<(zcat <female_reads1.fq>) <(zcat <female_reads2.fq>)

jellyfish histo -o male_reads.histo male_reads.counts
jellyfish histo -o fema_reads.histo fema_reads.counts
```
Notes: The two histo files can be then uploaded onto GenomeScope website (http://qb.cshl.edu/genomescope/) to estimate the genome size separately for female and male plants 


### Examine Genome Quality 

#### Remove potential contaminants from the initial assembly
```
sh qc_scripts/contaminants_dbbuild.sh -d <working dir> -n 16   

sh qc_scripts/contaminants_remove.sh -d <working dir> -g <genome file> -r solanum \
-R human,bacteria,viral -n 16 -i 80 -c 50
```
Notes: After first building the contaminant database, you need to manually edit "DeconSeqConfig.pm" (see the example file in 'qc_scripts') in deconseq2 package directory before running 'contaminants_remove.sh'

#### Examine te base call quality using female/male illumina reads by Referee 
```
sh qc_scripts/basecall_qc.sh -d <working dir> -g <genome file> \
-f <female_reads1.fq> -r <female_reads2.fq>

sh qc_scripts/basecall_qc.sh -d <working dir> -g <genome file> \
-f <male_reads1.fq> -r <male_reads2.fq>
```

#### Examine gene content using BUSCO against plant database
```
<path_to_BUSCO>/run_BUSCO.py -c 4 -i <genome file> -m genome -f \
-l <path_to_BUSCO>/embryophyta_odb9 -o busco_check
```


### Genome Annotation

#### Prepare species-specific repeats library
```
sh annotation_scripts/repeats_annotation.sh -d <working dir> -g <genome file> \
-T <Tpases020812DNA> -P <alluniRefprexp070416> -n 16 -S 1234
```

#### Prepare transcriptome fasta from RNA-seq reads (using HISAT2 and StringTie)
```
sh annotation_scripts/rnaseq_alignment.sh
```

#### Annotate genome using MAKER2 pipeline
```
maker -EXE   # generate maker_exe.ctl; need to add the PATH of required programs

sh annotation_scripts/gene_annotation.sh -d <working dir> -g <genome file> \
-e <transcriptome fasta> -p <uniprot solanaceae fasta> \
-R <Repeats.lib> -B <busco/embryophyta_odb9> \
-M maker_exe.ctl -I solanum_appendiculatum 
```

#### Check annotation statistics and pull out coding sequence with AED score <0.5
```
python annotation_scripts/annotation_by_aed.py --genome <genome file> --gff <GFF output from MAKER2> --AED 0.5
```

#### Add functional annotations to the predicted proteins
```
sh annotation_scripts/function_blast.sh -d <working dir> \
-s <proteins fasta from MAKER> -B <the directory saving all the blast datasets>

java -Xmx2g -jar <AHRD PATH>/ahrd.jar <your_use_case_file.yml>
```
Notes: Here the blast databasets include 'ITAG3.2_proteins.fasta', 'TAIR10_pep_20101214', and 'uniprot_sprot.fasta'. &nbsp; You need to modify the YAML file to add the PATH of BLAST results (see example file in 'annotation_scripts').

#### Infer putative syntenic blocks using Satsuma
```
SatsumaSynteny -t <tomato genome> -q <genome assembly> -o <output dir> -n 8
```


### Gene Family Analysis

#### Generate cluster input by using OrthoFinder
```
orthofinder -f <protein sequences directory>
```
Notes: the directory contains input protein fasta files, with one file per species.  For S. appendiculum, the fasta file was the one generated from 'annotation_by_aed.py'. &nbsp; The 'Orthogroups.GeneCount.csv' file in the OrthoFinder outputs can be then modified in the format a little bit and renamed as 'unfiltered_cafe_input.txt' (see the example file in 'cafe_scripts') 

#### run CAFE to infer gene family dynamics
```
python cafe_scripts/cafetutorial_clade_and_size_filter.py \
-i unfiltered_cafe_input.txt -o filtered_cafe_input.txt -s

sh cafe.sh
```
Notes: To set up CAFE, check more detail at https://iu.app.box.com/v/cafetutorial-pdf


### Transcriptome Analysis

#### Generate read count table by using FeatureCount 
```
featurecount -a <GTF file> -t exon -p -o raw_featurecount_data.txt <bam_file1> <bam_file2> <...>
```
Notes: The GFF file is generated above from 'annotation_by_aed.py'; &nbsp; Those BAM files are generated from 'rnaseq_alignment.sh'

#### Perform DEG analysis
```
Rscript DEG_scripts/DEG-analysis.R
```


### Sex Determination Analysis

#### Identify sex-specific kmers (SSKs)
```
sh sexy_scripts/sexy_kmers.sh
```

#### Identify illumina reads containing SSKs, aligned them to the genome, and count read depth across 10kb-windows
```
sh sexy_scripts/sexy_illumina.sh
```

#### Hijack the assembly script to generate corrected PacBio reads by sex
```
sh assembly_script/fema_assembly.sh
sh assembly_script/male_assembly.sh
```
Note: The goal is to obtain corrected pacbio long reads only. &nbsp; You can abort the script after you can find the reads fasta file "mr.41.15.17.0.029.1.fa" in the working directory.

#### Identify corrected pacbio reads (from separate genome assembly for the female/male) containing SSKs and aligned them to the genome
```
sh sexy_scripts/sexy_pacbio.sh
```

#### Examine population differentiation (Fst and Dxy) between two sexes
```
sh sexy_scripts/sexy_variants.sh
```

</br>
