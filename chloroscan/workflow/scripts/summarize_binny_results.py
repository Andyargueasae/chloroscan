import pandas as pd 
from os import listdir
import numpy as np
import sys
import logging
from Bio import SeqIO 
import urllib.request
import os
from pathlib import Path
from summarization_utils import *


# set up global variables that could work as input files.
NOT_BINNED = ""
depth_set = 0.0
CONTIG_ID = "contig id"
CONTIG_SEQ = "contig sequence"
CONTIG_DEP = "contig depth"
CONTIG_LEN = "contig length"
CONTIG_GC = "GC contents"
CONTIG_BIN = "Contig2Bin"
CONTIG_TAXON = "Taxon per Contig"
CONTIG_marker = "markers on the contig"
assemblyfasta = snakemake.params['assembly_files']
assemblydepth = snakemake.params['assembly_depth']
assembly_annot = snakemake.params['marker_gene_gff']
assembly_bin_dir = snakemake.params['bins']
taxonomy_names = snakemake.params['names_dump']
MMA_SUMMARY=snakemake.params['MMA_summary']
contig_level_annotation = os.path.join(snakemake.input['contig_level_annotation'], snakemake.params['ANNOT_FILE'])

# Here we want to know if the results from binny are empty.
if ((not Path(contig_level_annotation).exists() or os.stat(contig_level_annotation).st_size == 0) 
    and (not Path(assemblyfasta).exists()) 
    and (not Path(assemblydepth).exists()) 
    and (not Path(assembly_annot).exists())) :
    Path(contig_level_annotation).write_text("")
    print("No files supporting this rule to run were found, exit.\n")
    Path(snakemake.output[0]).write_text("")
    sys.exit(0)

# Meanwhile, if all those files are not there, we should simply go out.

dataset_bin_df = pd.DataFrame()
assembly_fasta_contigs = list(SeqIO.parse(assemblyfasta, "fasta"))
contig_id = []
contig_seq = []
contig_length = []
contig_gc = []
store_basic_info(contig_id, contig_seq, contig_length, contig_gc, assembly_fasta_contigs)
# formulate the length of the array for the following storage.
num_contigs = len(contig_id)

bin_array = list(np.repeat("", num_contigs))
bin_depth = list(np.repeat(0, num_contigs))
marker_per_contig = list(np.repeat("", num_contigs))
taxon_array = list(np.repeat("", num_contigs))

bin_files = listdir(assembly_bin_dir)

form_bin_array(assembly_bin_dir, bin_files, bin_array, contig_id)
form_marker_array(marker_per_contig, assembly_annot, contig_id)
form_depth_array(bin_depth, contig_id, assemblydepth)

form_contig_taxa_array(taxon_array, contig_id, taxonomy_names, contig_level_annotation)

dataset_bin_df[CONTIG_ID] = contig_id
dataset_bin_df[CONTIG_GC] = contig_gc
dataset_bin_df[CONTIG_DEP] = bin_depth
dataset_bin_df[CONTIG_LEN] = contig_length
dataset_bin_df[CONTIG_TAXON] = taxon_array
dataset_bin_df[CONTIG_marker] = marker_per_contig
dataset_bin_df[CONTIG_BIN] = bin_array
dataset_bin_df[CONTIG_SEQ] = contig_seq

# dataset_bin_df.to_csv(snakemake.output[0])

dataset_bin_df.to_csv(snakemake.output[0], sep="\t")

# Here we want to add some summary to MAGs.
bin_set=set(bin_array)
medium_bins=0
high_bins=0
total_bins=len(bin_set)-1
for i in bin_set:
    if i != "":
        file_name=i.replace(".fasta","")
        essential_info=file_name.split('_')
        # breakpoint()
        # Make sure you include both possibilities.
        if ("C" in essential_info[-3]) and ("P" in essential_info[-2]):            
            completeness=int(essential_info[-3].replace("C", ""))
            purity=int(essential_info[-2].replace("P", ""))
        #else, they are in fourth and third position.
        else:
            completeness=int(essential_info[-6].replace("C", ""))
            purity=int(essential_info[-5].replace("P", ""))

        if (completeness > 70) and (completeness < 90):
            medium_bins+=1
        else:
            high_bins+=1
    else:
        continue

with open(MMA_SUMMARY, "a") as MS:
    MS.write("Total number of MAGs: {}\n".format(total_bins))
    MS.write("Total number of medium-quality (completeness: 70% to 90%): {}\n".format(medium_bins))
    MS.write("Total number of high-quality (completeness: > 90%): {}\n".format(high_bins))
    

