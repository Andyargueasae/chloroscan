import pandas as pd 
from os import listdir
import numpy as np
import sys
import logging
from Bio import SeqIO 
import urllib.request
import os
from pathlib import Path


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
    # You wrote it empty. Then test it to be empty.
    # Meanwhile, the cross-ref.tsv should also be written.
    Path(snakemake.output[0]).write_text("")
    sys.exit(0)

# Meanwhile, if all those files are not there, we should simply go out.


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
        bin_path = str(bin_dir) + "/" + i 
        individual_bin = list(SeqIO.parse(bin_path, "fasta"))
        for elem in individual_bin:
            index_in_list = contig_id_list.index(elem.id)
            bin_name = i.replace(".fasta","")
            # print(bin_name)
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
        # print(key)
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
    # Here we missed the point that last year I only took 500 as length cutoff, if I test using other length cutoffs < 500, problems arise.
    for i in range(len(depth_lines)):
        contig, depth = depth_lines[i].replace("\n", "").split("\t",1)[0], depth_lines[i].replace("\n","").split("\t",1)[1]
        if contig in contig_id:
            empty_depth_array[contig_id.index(contig)] = depth
    return empty_depth_array

def store_basic_info(contig_id_array, contig_seq_array, contig_len_array, contig_gc_array, assembly_contigs):
    # simply change all of the files
    for i in assembly_contigs:
        contig_id_array.append(i.id)
        contig_seq_array.append(i.seq)
        contig_len_array.append(len(i.seq))
        contig_gc_array.append(GC_content(str(i.seq)))

    return
def find_finest_taxon(id_hierarchy, all_organisms):
    i = 0
    length_id_hierarchy = len(id_hierarchy)
    while i < length_id_hierarchy:
        try:
            that_organism = all_organisms[id_hierarchy[i]]
            i+=1
        except KeyError:
            break
    return all_organisms[id_hierarchy[i-1]]

def all_organisms_dict(names_dump):
    all_organisms = dict()
    with open(names_dump, "r") as nd:
        org_lines = nd.readlines()
    
    for i in org_lines:
        i = i.replace("\n", "").replace("\t", "")
        org_line = i.split("|")
        if "scientific name" in org_line:
            all_organisms[org_line[0]] = org_line[1]
    return all_organisms

def form_contig_taxa_array(empty_taxon_array, contig_id, names_dump, contig_annotation):
    # follow the same order as the contig_id.
    with open(contig_annotation, "r") as contig_f:
        contig_taxonomy_lines = contig_f.readlines()

    info_contig_taxonomy = contig_taxonomy_lines[1:]
    contig2taxid=dict()
    for i in info_contig_taxonomy:
        line_info = i.split("\t")
        if line_info[1] == "no taxid assigned":
            contig2taxid[line_info[0]] = "0"
        else:
            contig2taxid[line_info[0]] = line_info[-2]
    
    all_organisms = all_organisms_dict(names_dump)

    for i in contig2taxid.keys():
        index_contig = contig_id.index(i)
        if "*" in contig2taxid[i]:
            taxon_i = contig2taxid[i].replace("*","")
        taxon_i = contig2taxid[i]
        if taxon_i == "0":
            empty_taxon_array[index_contig] = "not classified"
        else:
            taxon_hierarchy = taxon_i.split(";")
            finest_classification = find_finest_taxon(taxon_hierarchy, all_organisms)
            empty_taxon_array[index_contig] = finest_classification

    return

dataset_bin_df = pd.DataFrame()
# breakpoint()
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
    

