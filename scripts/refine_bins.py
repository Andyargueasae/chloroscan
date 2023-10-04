import Bio
from Bio import SeqIO
import os
from os import listdir
import pandas as pd
import numpy as np

# Some global variable inputs.
output_dir_for_refined_bins = snakemake.output[0]

# Inputs and params.
cross_reference_table = snakemake.input['cross_ref']
original_bins = snakemake.input['original_bins']+"/bins"
CAT_prediction = snakemake.input['CAT_prediction'] + "/out.CAT.contig2classification.txt"
BACTERIA_ID = snakemake.params['BACTERIA_ID']

#read in some files.
os.mkdir(output_dir_for_refined_bins)

# breakpoint()
if (os.stat(cross_reference_table).st_size == 0):
    print("Empty cross reference table, exit.\n")
    sys.exit(0)


cross_ref_refine_df = pd.read_csv(cross_reference_table,sep="\t", index_col=0)
bin_list = listdir(original_bins)
orig_CAT_df = pd.read_csv(CAT_prediction, sep="\t")
lineages = orig_CAT_df['lineage']

list_of_seqrec = dict()
for i in bin_list:
    filename = original_bins + "/" + i
    list_of_seqrec[filename] = list(SeqIO.parse(filename, "fasta"))


# First, identify which contig is bacterial and has no marker genes on it.
refined_contigs_dict = dict()

to_be_modified_id = []
for i in list_of_seqrec.keys():
    contig_comp = list_of_seqrec[i]
    contig_ids = [i.id for i in contig_comp]
    print("contig composition before refinement:\n", contig_ids)
    print("length equals to {}".format(len(contig_ids)))
    for elem in contig_comp:
#         print(elem.id)
        elem_id = elem.id
        elem_in_crossref = cross_ref_refine_df.loc[cross_ref_refine_df['contig id'] == elem_id]
#         print(elem_in_crossref)
        elem_row = orig_CAT_df.loc[orig_CAT_df["# contig"] == elem.id]
        elem_lineage = np.array(elem_row['lineage'])[0]
#         print(elem_in_crossref)
        bacteria_taxonomy_id = BACTERIA_ID
        if not pd.isna(elem_lineage):
            elem_lineage_hierarchy = elem_lineage.split(";")
            print(elem_lineage_hierarchy)
            print(np.array(elem_in_crossref['markers on the contig'])[0])
            if (bacteria_taxonomy_id in elem_lineage_hierarchy) and (pd.isna(np.array(elem_in_crossref['markers on the contig'])[0])):
                print("This is a bacterial contamination.")
                contig_ids.remove(elem_id)
                to_be_modified_id.append(elem_id)
            else:
                continue
        else:
            continue
    print("contig composition after refinement: \n", [j for j in contig_ids])
    print("length equals to {}".format(len([l for l in contig_ids])))
    refined_contig_ids = contig_ids
    
    # Now store new bag of contigs.
    refined_contigs_dict[i] = []
    for req in contig_comp:
        if req.id in refined_contig_ids:
            refined_contigs_dict[i].append(req)
    print("Number of contigs from the refined bin {}".format(len(refined_contigs_dict[i])))
    print("--------------------------------------------------")

# update the df, and change it back to the place!
for id_change in to_be_modified_id:
    row = cross_ref_refine_df.loc[cross_ref_refine_df['contig id'] == id_change]
    #get the index.
    row_index = row.index[0]
    cross_ref_refine_df.at[row_index, "Contig2Bin"] = ""

# refined_contigs_dict
for i in refined_contigs_dict.items():
    filename = output_dir_for_refined_bins + "/" + "refined_" + i[0].split("/")[-1]
    print(filename)
    with open(filename, "w") as refined:
        for j in i[-1]:
            print(j.id)
            refined.write(">")
            refined.write(j.id)
            refined.write("\n")
            refined.write(str(j.seq))
            refined.write("\n")


cross_ref_refine_df.to_csv(snakemake.input['cross_ref'], sep="\t")
del cross_ref_refine_df