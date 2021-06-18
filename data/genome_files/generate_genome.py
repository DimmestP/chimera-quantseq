# -*- coding: utf-8 -*-
"""
Created on Tue May 18 00:59:06 2021

@author: Weronika
"""


import pandas as pd
from math import isnan
import subprocess

GTF_COLNAMES = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
FEATURES = ["Promoter", "CDS", "Terminator", "URA3", "M1", "M2", "M3"]

def loadGtfFile(file_path):
    gtf_df = pd.read_csv(file_path, sep="\t", header=None, comment="#")
    gtf_df.columns = GTF_COLNAMES
    return gtf_df

def extractGeneIdFromGtfAttributes(gtf_df):
    gtf_df[['gene_id','attributes']] = gtf_df['attributes'].str.split(';', 1, expand=True)
   
    # remove "gene_id " at the beginning of the field
    gtf_df["gene_id"] = gtf_df["gene_id"].str.strip("gene_id ")
    gtf_df["gene_id"] = gtf_df["gene_id"].str.strip("\"")
    return gtf_df

def findFeatureCoordinatesInGtf(gene_id, feature, gtf_df):
    gtf_df_sub = gtf_df[(gtf_df["gene_id"] == gene_id) & (gtf_df["feature"] == feature)]
    
    chromosome = gtf_df_sub["seqname"].iloc[0]
    start = gtf_df_sub["start"].iloc[0]
    stop = gtf_df_sub["end"].iloc[0]
    strand = gtf_df_sub["strand"].iloc[0]

    return chromosome, start, stop, strand


# generate annotations for 3'-UTRs of reference genes
# load csv file with ref genes and max 3' UTR lengths
reference_genes_utr = pd.read_csv("./chosen_terminator_max_3UTR_length.csv", sep=",")
print(reference_genes_utr.head(5))

# load Scer gtf file to extract reference genes coordinates
annotation_gtf = loadGtfFile("./genomes/GCF_000146045.2_R64_genomic.gtf")
annotation_gtf_plus_gene = extractGeneIdFromGtfAttributes(annotation_gtf)
print(annotation_gtf_plus_gene.head(5))

ref_genes = reference_genes_utr["transcript_name"].tolist()
ref_utr_lengths = reference_genes_utr["UTR3_length"]
ref_gene_annotations = []
for ref_gene, ref_utr_length in zip(ref_genes, ref_utr_lengths):
    gtf_per_gene = annotation_gtf_plus_gene[annotation_gtf_plus_gene["gene_id"] == ref_gene]
    print(gtf_per_gene)
    chromosome, start, stop, strand = findFeatureCoordinatesInGtf(ref_gene, "stop_codon", gtf_per_gene)
    print(chromosome, start, stop, strand)
    
    # find coordinates of 3'-UTR
    if strand == "+":
        utr_start = stop + 1
        utr_stop = stop  + ref_utr_length
    else:
        utr_stop = start - 1
        utr_start = start - ref_utr_length
                    
    attributes_str = "gene_id \"" + ref_gene + "\";"

    new_annotation_str = "\t".join([chromosome, "chimera-project", "utr-3p", str(int(utr_start)), str(int(utr_stop)), ".", strand, ".", attributes_str]) 
    
    ref_gene_annotations.append(new_annotation_str)
    
print(ref_gene_annotations)


# generate plasmid features annotations (prom, CDS, term, motif inserts, marker)
plasmid_features = pd.read_csv("./plasmid_features.csv", sep = ";")
print(plasmid_features.head(5))


construct_names = plasmid_features['Construct_name'].tolist()

for construct_name in construct_names:
    with open('./genomes/GCF_000146045.2_R64_genomic.gtf') as fp:
        data = fp.read()

    for ref_gene_annotation in ref_gene_annotations:
        data += ref_gene_annotation
        data += "\n"
        
    for feature in FEATURES:
        feat_start_col = feature + "_start"
        feat_stop_col = feature + "_stop"
        
        feat_start = plasmid_features[plasmid_features["Construct_name"] == construct_name][feat_start_col].iloc[0]
        feat_stop = plasmid_features[plasmid_features["Construct_name"] == construct_name][feat_stop_col].iloc[0]
        if isnan(feat_start) == False:

            if feature == "URA3":
                strand = "-"
            else:
                strand = "+"
                
            attributes_str = "gene_id \"" + construct_name + "\";"
            annotation_str = "\t".join([construct_name, "chimera-project", feature, str(int(feat_start)), str(feat_stop), ".", strand, ".", attributes_str])
            data += annotation_str
            data += "\n"


    # sorting GTF file - requires gff3sort package
    with open ('./genomes/GCF_000146045.2_R64_genomic_' + construct_name + '_unsorted.gtf', 'w') as fp:
        fp.write(data)
        
    command = ["gff3sort.pl",
               "--precise",
               "--chr_order",
               "natural",
               './genomes/GCF_000146045.2_R64_genomic_' + construct_name + '_unsorted.gtf']
    
    with open('./genomes/GCF_000146045.2_R64_genomic_' + construct_name + '.gtf', "w") as outfile:
        subprocess.run(command, stdout=outfile)
        

# add plasmid maps to genome

for construct_name in construct_names:
    fasta_file = "./linearized_plasmids/linear_" + construct_name +".fa"
    genome_file = "./genomes/GCF_000146045.2_R64_genomic_" + construct_name +".fa"

    command = ["cat",
               "./genomes/GCF_000146045.2_R64_genomic.fna",
               fasta_file]
    
    with open(genome_file, "w") as outfile:
        subprocess.run(command, stdout=outfile)
        
# generate IGV genome files (which are zip files with ".genome" extension)
# NOTE: If .genome file already exists, new files will be added or replaced, but old files might remain - clean up the dir before running
    
for construct_name in construct_names:
    genome_file = "./genomes/GCF_000146045.2_R64_genomic_" + construct_name +".fa"
    gtf_file = './genomes/GCF_000146045.2_R64_genomic_' + construct_name + '.gtf'
    
    command = ["samtools", "faidx",
               genome_file]
    subprocess.run(command)
    
    property_file = ("fasta=true\n" \
                    "fastaDirectory=false\n" \
                    "ordered=true\n" \
                    "id=" + construct_name + "_genome\m" \
                    "name=" + construct_name + "_genome\n" \
                    "geneFile=" + "GCF_000146045.2_R64_genomic_" + construct_name + ".gtf\n" \
                    "sequenceLocation=" + "GCF_000146045.2_R64_genomic_" + construct_name +".fa\n")
 
    with open ('./genomes/property.txt', 'w') as fp:
        fp.write(property_file)
        
    command = ["zip", "-j", "-r",
               "./genomes/" + construct_name + ".genome",
               "./genomes/property.txt",
               gtf_file]
    subprocess.run(command)   
    
    command = ["rm", "./genomes/property.txt"]
    subprocess.run(command)