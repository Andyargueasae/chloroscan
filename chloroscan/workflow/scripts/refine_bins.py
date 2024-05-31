import Bio
from Bio import SeqIO
import os
from os import listdir
import pandas as pd
import numpy as np
import sys
# Do we need to apply another constraint: plastid bin lengths?

# Some global variable inputs. This script needs to be better ellaborated.
output_dir_for_refined_bins = snakemake.output[0]

# Inputs and params.
cross_reference_table = snakemake.input['cross_ref']
original_bins = snakemake.input['original_bins']+"/bins"
CAT_prediction = snakemake.input['CAT_prediction'] + "/out.CAT.contig2classification.txt"
BACTERIA_ID = snakemake.params['BACTERIA_ID']

#read in some files.
os.mkdir(output_dir_for_refined_bins)

## breakpoint()
if (os.stat(cross_reference_table).st_size == 0):
    print("Empty cross reference table, exit.\n")
    sys.exit(0)

# If corgi has the outputs while binny doesn't, we stop the refinement.
# if os.path.exists(original_bins):
#     binny_check = listdir(original_bins)
#     if len(binny_check) == 0:
#         print("No bins produced, system shut down.")
#         sys.exit(0)
#     else:
#         print("Ok, your bins are good to go.")


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
        elem_id = elem.id
        elem_in_crossref = cross_ref_refine_df.loc[cross_ref_refine_df['contig id'] == elem_id]
        elem_row = orig_CAT_df.loc[orig_CAT_df["# contig"] == elem.id]
        elem_lineage = np.array(elem_row['lineage'])[0]
        bacteria_taxonomy_id = BACTERIA_ID
        if not pd.isna(elem_lineage):
            elem_lineage_hierarchy = elem_lineage.split(";")
            print(elem_lineage_hierarchy)
            print(np.array(elem_in_crossref['markers on the contig'])[0])
            if (bacteria_taxonomy_id in elem_lineage_hierarchy) and (pd.isna(np.array(elem_in_crossref['markers on the contig'])[0])):
                print("This is a bacterial contamination.")
                contig_ids.remove(elem_id)
                to_be_modified_id.append(elem_id)
            # If contig is not even identified as bacteria and you did not have any markers, remove.
            # But the only problem is the off-target: if the contig has some undocumented markers and was not 
            # classified algal sequence by CAT, it will be removed as well.

            # We will have to go back to the
            elif (pd.isna(np.array(elem_in_crossref['markers on the contig'])[0])):
                print("This is a contamination.")
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

def refine_bins_renaming(seqs, cross_ref_df):
    taxonomy_dict = dict()
    # breakpoint()
    seq_ids = [i.id for i in seqs]
    print(seq_ids)
    for seq_id in seq_ids:
        the_row = cross_ref_df.loc[cross_ref_df['contig id'] == seq_id]
        the_taxon = the_row['Taxon per Contig'].item()
        the_length = the_row['contig length'].item()
        if the_taxon not in taxonomy_dict.keys():
            taxonomy_dict[the_taxon] = 0
        taxonomy_dict[the_taxon] += the_length
    
    print(taxonomy_dict)
    sorted_dict_items = sorted(taxonomy_dict.items(), key = lambda x:x[1])
    print(sorted_dict_items)
    taxon_with_longest_length = sorted_dict_items[-1][0]
    return taxon_with_longest_length
    
for i in refined_contigs_dict.items():
    filename = output_dir_for_refined_bins + "/" + "refined_" + i[0].split("/")[-1]    
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