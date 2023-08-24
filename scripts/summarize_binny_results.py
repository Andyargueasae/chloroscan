import pandas as pd 
from os import listdir
import numpy as np
import sys
import logging
from Bio import SeqIO 

# set up global variables that could work as input files.
# You may also work out the visualization of the individual dataset info here. 
# Note: this time, the individual contig identity shall not be put into the process.

NOT_BINNED = ""
depth_set = 0.0
CONTIG_ID = "contig id"
CONTIG_SEQ = "contig sequence"
CONTIG_DEP = "contig depth"
CONTIG_LEN = "contig length"
CONTIG_GC = "GC contents"
CONTIG_BIN = "Contig2Bin"
CONTIG_marker = "markers on the contig"
#Assign global variables from snakemake.
# @click.command()
# @click.argument("bins_directory", type = click.Path(exists=True))
# @click.option("--assemblyfasta", type = click.Path(exists=True))
# @click.option("--assemblydepth", type = click.Path(exists=True))
# @click.option("--assemblygff", type = click.Path(exists=True))
assemblyfasta = snakemake.params['assembly_files']
assemblydepth = snakemake.params['assembly_depth']
assembly_annot = snakemake.params['marker_gene_gff']
assembly_bin_dir = snakemake.input[0]


def GC_content(sequence):
    G_count = sequence.count('G')
    C_count = sequence.count('C')
    return (G_count + C_count)/len(sequence)

def form_bin_array(bin_dir, collected_bins, bin_array, contig_id_list):
    # bin_dir: bin_directory inside binny;
    # collected_bins: bin_file name;
    # bin_array: the empty list containing len(contig_id) elements.
    # contig_id_list: the list of contigs with their id.
    
    for i in collected_bins:
        bin_path = bin_dir + "/" + i 
        individual_bin = list(SeqIO.parse(bin_path, "fasta"))
        for elem in individual_bin:
            index_in_list = contig_id_list.index(elem.id)
            bin_name = i.replace(".fasta","")
            print(bin_name)
            bin_array[index_in_list] = bin_name
    return 

def form_marker_array(empty_marker_array, annot_file, contig_id):
    # empty_marker_Array needs to be defined before.
    with open(annot_file) as af:
        annot_lines = af.readlines()
    marker_dict = {}
    default = "no markers in the database"
    # len(annot_lines)
    for i in range(len(annot_lines)):
        i_th_ORF = annot_lines[i].split("\t")
    #     print(i_th_ORF)
        
        key = i_th_ORF[0]
        ORF = i_th_ORF[-1]
        print(key)
    #     print(ORF)
        ORF_split = ORF.replace("\n","").split(";")
        # print(ORF_split[-1])
        # Initiate the gathering.
        if ("checkm_marker" in ORF_split[-1]) and (key not in marker_dict.keys()):
            marker_dict[key] = ""
            marker_dict[key] += ORF_split[-1].split("=")[-1]
        elif ("checkm_marker" in ORF_split[-1]) and (key in marker_dict.keys()):
            marker_dict[key] += ","
            marker_dict[key] += ORF_split[-1].split("=")[-1]
        else:
            if key not in marker_dict.keys():
                marker_dict[key] = ""
            else:
                continue
    
    for i in marker_dict.items():
        key = i[0]
        value = i[-1]
        index = contig_id.index(key)
        empty_marker_array[index] = value

    # Now the marker array is no longer empty.
    return empty_marker_array

def form_depth_array(empty_depth_array, contig_id, depth_file):
    with open(depth_file) as df:
        depth_lines = df.readlines()
    
    for i in range(len(depth_lines)):
        contig, depth = depth_lines[i].split("\t")[0], depth_lines[i].split("\t")[1]
        empty_depth_array[i] = depth
    return empty_depth_array

def store_basic_info(contig_id_array, contig_seq_array, contig_len_array, contig_gc_array, assembly_contigs):
    # simply change all of the files
    for i in assembly_contigs:
        contig_id_array.append(i.id)
        contig_seq_array.append(i.seq)
        contig_len_array.append(len(i.seq))
        contig_gc_array.append(GC_content(str(i.seq)))

    return

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

bin_files = listdir(assembly_bin_dir)

form_bin_array(assembly_bin_dir, bin_files, bin_array, contig_id)
form_marker_array(marker_per_contig, assembly_annot, contig_id)
form_depth_array(bin_depth, contig_id, assemblydepth)


dataset_bin_df[CONTIG_ID] = contig_id
dataset_bin_df[CONTIG_GC] = contig_gc
dataset_bin_df[CONTIG_DEP] = bin_depth
dataset_bin_df[CONTIG_LEN] = contig_length
dataset_bin_df[CONTIG_marker] = marker_per_contig
dataset_bin_df[CONTIG_BIN] = bin_array
dataset_bin_df[CONTIG_SEQ] = contig_seq


dataset_bin_df.to_excel(snakemake.output[0])
