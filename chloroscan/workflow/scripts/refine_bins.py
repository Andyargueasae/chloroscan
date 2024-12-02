import Bio
from Bio import SeqIO
import os
from os import listdir
import pandas as pd
import numpy as np
import sys
import logging

def refine_bins_renaming(seqs, cross_ref_df):
    taxonomy_dict = dict()
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

output_dir_for_refined_bins = snakemake.output[0]
output_summary_info = snakemake.output[1]

# Inputs and params.

cross_reference_table = snakemake.input['cross_ref']
original_bins = snakemake.input['original_bins']+"/bins"
CAT_prediction = snakemake.input['CAT_prediction'] + "/out.CAT.contig2classification.txt"
BACTERIA_ID = snakemake.params['BACTERIA_ID']
EUKARYA_ID = "2759"
ROOT_ID = "1"
log_file = snakemake.log[0]
logging.basicConfig(
    filename=log_file,  # Specify the log file
    filemode='a',                # 'a' for append mode, 'w' for write (overwrite)
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Log format
    level=logging.INFO           # Set the minimum log level to INFO
)

logger = logging.getLogger(__name__)

os.mkdir(output_dir_for_refined_bins)

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
    logging.info(f"Inspect bin {i.split('/')[-1]} contig composition.")
    contig_comp = list_of_seqrec[i]
    contig_ids = [i.id for i in contig_comp]
    contig_length = [len(i.seq) for i in contig_comp]
    logging.info(f"contig composition before refinement: {','.join(contig_ids)}\n")
    logging.info(f"Number of contigs: {len(contig_ids)}")
    logging.info(f"Summed Length of MAG: {sum(contig_length)}")
    # Length of MAG should be reported.
    for elem in contig_comp:
        elem_id = elem.id
        elem_in_crossref = cross_ref_refine_df.loc[cross_ref_refine_df['contig id'] == elem_id]
        elem_row = orig_CAT_df.loc[orig_CAT_df["# contig"] == elem.id]
        elem_lineage = np.array(elem_row['lineage'])[0]
        bacteria_taxonomy_id = BACTERIA_ID
        if not pd.isna(elem_lineage):
            elem_lineage_hierarchy = elem_lineage.split(";")
            logging.info("CAT/BAT taxonomic assignment of contig {}: {}".format(elem_id, elem_lineage_hierarchy))
            logging.info("Marker genes found on the contig: {}".format(np.array(elem_in_crossref['markers on the contig'])[0]))
            if (bacteria_taxonomy_id in elem_lineage_hierarchy) and (pd.isna(np.array(elem_in_crossref['markers on the contig'])[0])):
                logging.info("contig id: {} This is a bacterial contamination.".format(elem_id))
                contig_ids.remove(elem_id)
                to_be_modified_id.append(elem_id)
                with open(output_summary_info, "a") as si:
                    si.write(f"MAG {i.split('/')[-1]} has suspicious contamination contig: {elem_id}. Taxonomy: {elem_lineage_hierarchy}.")

            elif (EUKARYA_ID in elem_lineage_hierarchy) and (pd.isna(np.array(elem_in_crossref['markers on the contig'])[0])):
                logging.info("Contig id: {}. This is a eukaryotic contig without marker genes on it.".format(elem_id))
                contig_ids.remove(elem_id)
                to_be_modified_id.append(elem_id)
                with open(output_summary_info, "a") as si:
                    si.write(f"MAG {i.split('/')[-1]} has a eukaryotic contig without single-copy marker genes in database: {elem_id}. Taxonomy: {elem_lineage_hierarchy}.")

            elif (ROOT_ID in elem_lineage_hierarchy) and (pd.isna(np.array(elem_in_crossref['markers on the contig'])[0])):
                logging.info("Contig id: {}. This is a contig with ambiguous taxonomic classification and without recognizable A2K marker genes, discarded.".format(elem_id))
                contig_ids.remove(elem_id)
                to_be_modified_id.append(elem_id)
                with open(output_summary_info, "a") as si:
                    si.write(f"MAG {i.split('/')[-1]} has an unclassified contig: {elem_id}.")

            else:
                continue
        else:
            continue
    logging.info(f"contig composition after refinement: {','.join(contig_ids)}\n")
    logging.info(f"length equals to {len(contig_ids)}")
    refined_contig_ids = contig_ids
    
    # Now store new bag of contigs.
    refined_contigs_dict[i] = []
    for req in contig_comp:
        if req.id in refined_contig_ids:
            refined_contigs_dict[i].append(req)
    logging.info("Number of contigs from the refined bin {}".format(len(refined_contigs_dict[i])))
    # Length of MAG should be reported.
    new_contig_length = [len(i.seq) for i in contig_comp if i.id in contig_ids]
    logging.info(f"MAG length after refinement: {sum(new_contig_length)}")
    logging.info("-"*50)

# Do we want to update the df?
for id_change in to_be_modified_id:
    row = cross_ref_refine_df.loc[cross_ref_refine_df['contig id'] == id_change]
    row_index = row.index[0]
    cross_ref_refine_df.at[row_index, "Contig2Bin"] = ""
    
for i in refined_contigs_dict.items():
    filename = output_dir_for_refined_bins + "/" + "refined_" + i[0].split("/")[-1]    
    with open(filename, "w") as refined:
        for j in i[-1]:
            logging.info(j.id)
            refined.write(">")
            refined.write(j.id)
            refined.write("\n")
            refined.write(str(j.seq))
            refined.write("\n")


cross_ref_refine_df.to_csv(snakemake.input['cross_ref'], sep="\t")
