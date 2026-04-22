import pandas as pd 
from os import listdir
import numpy as np
import sys
from Bio import SeqIO 
import os
from pathlib import Path
from summarization_utils import *
import argparse


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

def main():
    parser= argparse.ArgumentParser(description="Summarize the binning results from binny and add taxonomic annotation to the contigs.")
    parser.add_argument("--assemblyfasta", help="The assembly fasta file used for binning.")
    parser.add_argument("--assemblydepth", help="The depth file of the assembly used for binning.")
    parser.add_argument("--assembly_annot", help="The marker gene annotation file of the assembly used for binning.")
    parser.add_argument("--assembly_bin_dir", help="The directory containing the bins from binny.")
    parser.add_argument("--MMA_SUMMARY", help="The output file for the summary of medium and high quality MAGs.")
    parser.add_argument("--contig_level_annotation", help="The contig level annotation file for the contigs.")
    parser.add_argument("--output", help="The output file for the summarized contig level annotation.")
    args = parser.parse_args()

    assemblyfasta = args.assemblyfasta
    assemblydepth = args.assemblydepth
    assembly_annot = args.assembly_annot
    assembly_bin_dir = args.assembly_bin_dir
    MMA_SUMMARY = args.MMA_SUMMARY
    contig_level_annotation = args.contig_level_annotation
    output = args.output

    if ((not Path(contig_level_annotation).exists() or os.stat(contig_level_annotation).st_size == 0) 
        and (not Path(assemblyfasta).exists()) 
        and (not Path(assemblydepth).exists()) 
        and (not Path(assembly_annot).exists())) :
        Path(contig_level_annotation).write_text("")
        print("No files supporting this rule to run were found, exit.\n")
        Path(output).write_text("")
        sys.exit(0)


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
    form_contig_taxa_array(taxon_array, contig_id, contig_level_annotation)

    dataset_bin_df[CONTIG_ID] = contig_id
    dataset_bin_df[CONTIG_GC] = contig_gc
    dataset_bin_df[CONTIG_DEP] = bin_depth
    dataset_bin_df[CONTIG_LEN] = contig_length
    dataset_bin_df[CONTIG_TAXON] = taxon_array
    dataset_bin_df[CONTIG_marker] = marker_per_contig
    dataset_bin_df[CONTIG_BIN] = bin_array
    dataset_bin_df[CONTIG_SEQ] = contig_seq

    dataset_bin_df.to_csv(output, sep="\t")

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
        
    return


if __name__ == "__main__":
    main()