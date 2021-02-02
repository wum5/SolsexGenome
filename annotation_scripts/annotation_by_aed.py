#!/usr/bin/env python
"""\n This program is to extract high-confidence gene annotions from MAKER2 \
output based on the annotation edit distance (AED). AED = 1 - AC, where \
AC = (SN + SP) / 2. SN = TP / (TP + FN); SP = TP / (TP + FP); The program \
gffread must be added into PATH"""


import os
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO



def extract_value(string):
    """
    Extract value after '='
    """
    value = string.split('=')[1]
    return value



def sort_geneid(old_ids, gene_pos, scaff_id):
    """
    Sort the gene ids by on on their coordinates
    """
    renamed_ids = {}
    gene_order = [i[0] for i in sorted(enumerate(gene_pos), key=lambda x:x[1])]
	
    for x in xrange(len(gene_order)):
        if x < 10:
            new_ids = scaff_id + 'g000' + str(x*10+10)
        elif 10 <= x < 100:
            new_ids = scaff_id + 'g00' + str(x*10+10)
        elif 100 <= x < 1000:
            new_ids = scaff_id + 'g0' + str(x*10+10)
        elif 1000 <= x < 10000:
            new_ids = scaff_id + 'g' + str(x*10+10)

        renamed_ids[old_ids[x]] = new_ids
	
    return renamed_ids



def rename_geneid(raw_gff3, updated_gff3):
    """
    Rename the genes based on their coordinates on scaffolds, and output generated is 
    the mapping between old ids and new ids
    """
    infile = open(raw_gff3, "r")
    outfile = open(updated_gff3, "w")
	
    geneids = {}
    prev_scaf = ''
    prev_record = ''
    gene_pos = []
    new_ids = []
    old_ids = []	
		
    for line in infile:
        line = line.rstrip()
        if '#' in line: continue
				
        if line.split('\t')[2] == 'gene':
            curr_scaf = line.split('\t')[0]
            gene_start = line.split('\t')[3]
            geneid = ((line.split('\t')[8]).split(';')[0]).split('ID=')[1]
            #print curr_scaf, prev_scaf
			
            if curr_scaf == prev_scaf:
                gene_pos.append(int(gene_start))
                old_ids.append(geneid)
                prev_record += line+'\n'

            else:
                geneids = sort_geneid(old_ids, gene_pos, prev_scaf)
                for key in geneids:
                    prev_record = prev_record.replace(key, geneids[key])
                prev_record = prev_record.replace('-mRNA-', '.')
                outfile.write(prev_record)
				
                prev_record = ''
                prev_scaf = curr_scaf
				
                gene_pos = []
                old_ids = []
                new_ids = []
                gene_pos.append(int(gene_start))
                old_ids.append(geneid)
                prev_record += line+'\n'
				
        else:
            prev_record += line+'\n'

    geneids = sort_geneid(old_ids, gene_pos, prev_scaf)
    for key in geneids:
        prev_record = prev_record.replace(key, geneids[key])
	
    prev_record = prev_record.replace('-mRNA-', '.')
    outfile.write(prev_record)
			
    infile.close()		
    outfile.close()	

	

def longest_isoform(infile):
    """
    remove redundant transcripts from a same gene based on the length of sequences, retaining
    each gene with the number of amino acids (AA) or bases (DNA) from the longest isoform
    """
    length_dict = {}
    index = infile.split('.')[0]
	
    for record in SeqIO.parse(infile, "fasta"):
        if 'mRNA' in record.id:
            name = record.id.split('-mRNA-')[0]
        elif '.' in record.id:
            name = record.id.split('.')[0]
			
        if name not in length_dict:
            length_dict[name]=len(str(record.seq))
        elif len(str(record.seq)) > length_dict[name]:
            length_dict[name]=len(str(record.seq))			
	
    outfile = open(index+"_NR.fasta", "w")
	
    for record in SeqIO.parse(infile, "fasta"):
        if 'mRNA' in record.id:
            name = record.id.split('-mRNA-')[0]
        elif '.' in record.id:
            name = record.id.split('.')[0]
			
        size = len(str(record.seq))
        if length_dict[name] == size:
            outfile.write(">"+record.id+'\n')
            outfile.write(str(record.seq)+'\n')
            length_dict[name] = -1
			
    outfile.close()
	
	
	
def gff3fasta(genome_seq, gff3file):
    """
    Extract CDS/protein sequences based on GFF3 file
    """
    index = (os.path.basename(gff3file)).split('.')[0]
    cmd="gffread -x %s_cds.fasta -g %s %s" % (index, genome_seq, gff3file)
    os.system(cmd)
    cmd="gffread -y %s_proteins.fasta -g %s %s" % (index, genome_seq, gff3file)
    os.system(cmd)
	
    longest_isoform(str(index)+'_cds.fasta')
    longest_isoform(str(index)+'_proteins.fasta')
	
	

def parse_maker_gff3(gff_file, AED_cutoff):
    """
    Extract high-confidence gene models from MAKER2 output using the AED cutoffs and 
    generate two output files: the updated GFF3 file and the non-redundant protein fasta
    """		
    index = gff_file.split('.')[0]
    infile = open(gff_file, "r")
	
    # initialization
    output = []
    outfile = open("tempfile", "w")
    gene_num, mRNA_num = 0, 0
    AED_thres={0.0:0, 0.1:0, 0.2:0, 0.3:0, 0.4:0, 0.5:0, 0.6:0, 0.7:0, 0.8:0, 0.9:0, 1.0:0}
	
    for line in infile:
        if 'maker' not in line: continue
        line = line.rstrip()
        feature = line.split('\t')[2]
        annotation = line.split('\t')[8]
		
        if feature == 'gene':
            # print the previous qualified gene model
            if len(output) > 1:
                gene_num += 1
                for element in output:
                    outfile.write(element+'\n')
					
            # re-initialization
            output = []
            output.append(line)
            mrna_ids = []
			
        elif feature == 'mRNA':
            AED = float(extract_value(annotation.split(';')[3]))
			
            if AED <= AED_cutoff:
                tmp = extract_value(annotation.split(';')[0])
                mrna_ids.append(tmp)
                mRNA_num += 1
                output.append(line)
					
            for key in AED_thres:
                if AED < key:
                    AED_thres[key]+=1
									
				
        elif feature in ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']:
            tmp = extract_value(annotation.split(';')[1])
            parents = tmp.split(",")
            for item in parents:
                if item in mrna_ids:
                    output.append(line)
									
        else: 
            pass
		
	
    # print the last qualified gene model
    if len(output) > 1:
        gene_num += 1
        for element in output:
            outfile.write(element+'\n')
	
    outfile.close()
    infile.close()
	
	
    # print the number genes matching the criteria to the screen		
    if AED_cutoff < 1:
        with open("genes_summary.txt", "w") as outfile:
            outfile.write("AED cutoff: %.2f \nNumber of genes: %d \nNumber of transcripts: %d\n" 
            % (AED_cutoff, gene_num, mRNA_num))

    plt.plot(AED_thres.keys(), AED_thres.values(), 'ro', label="AED")
    plt.xticks(AED_thres.keys())
    plt.ylabel('Number of Gene Models', fontsize=12)
    plt.xlabel('AED scores', fontsize=12)
    plt.title("Annotation quality check with MAKER's AED score", fontsize=15)
    plt.savefig("AED_plot.pdf")
			
	
    # write the updated GFF3 file
    updated_gff3 = index+'.aed.'+str(AED_cutoff)+'.gff3'	
    rename_geneid('tempfile', updated_gff3)
    os.system('rm tempfile')	

    return updated_gff3
	


def gff_stats(gff_file):
    """
    Obtain the statistics of gene annotations from the updated gff file
    """
    exon_num, gene_num, exon_size, cds_size, gene_size = 0, 0, 0, 0, 0

    with open(gff_file, "r") as infile:
        for line in infile:
            if 'maker' not in line: continue
            line = line.rstrip()
            feature = line.split('\t')[2]
            start = int(line.split('\t')[3])
            end = int(line.split('\t')[4])
            size = end - start + 1	
		
            if feature == 'mRNA':
                gene_num += 1
                gene_size += size
			
            elif feature == 'exon':
                exon_num += 1
                exon_size += size
		
            elif feature == 'CDS':
                cds_size += size
		
    avg_exon_num = float(exon_num)/gene_num
    avg_exon_size = float(exon_size)/exon_num
    avg_cds_size = float(cds_size)/gene_num
    avg_transcript_size = float(exon_size)/gene_num
    intron_size = float(gene_size-exon_size)/(exon_num-gene_num)

    with open("genes_summary.txt", "a") as outfile:
        outfile.write("Mean number of exons per gene: %f\n" % (avg_exon_num))
        outfile.write("Mean exon size: %f\n" % (avg_exon_size))
        outfile.write("Mean transcript length: %f\n" % (avg_transcript_size))
        outfile.write("Mean CDS size: %f\n" % (avg_cds_size))
        outfile.write("Mean intron size: %f\n" % (intron_size))
	
	
		

def main():
    parser = argparse.ArgumentParser(prog="annotation_by_AED.py", description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--genome', help=("input genome sequence in fasta"), required=True)
    parser.add_argument('--gff', help=("the raw GFF3 file from MAKER output"), required=True)
    parser.add_argument('--AED', help=("the cutoff AED score used for filtering"), default=1)
    args = parser.parse_args()
    
    gff_file = args.gff
    genome_seq = args.genome
    AED_cutoff = float(args.AED)

    updated_gff3 = parse_maker_gff3(gff_file, AED_cutoff)
    gff3fasta(genome_seq, updated_gff3)
    gff_stats(updated_gff3)
	
	

if __name__ == "__main__":
    main()
	
