from Bio import SeqIO
import os
from os import listdir
import pandas as pd
import numpy as np
import sys
import logging
from pathlib import Path
import argparse

BACTERIA_ID = "2"
EUKARYA_ID = "2759"
ROOT_ID = "1"

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

def main():
    parser = argparse.ArgumentParser(description="Refine the bins from binny by removing contigs that are likely to be contaminants based on their taxonomic annotation and the presence of marker genes.")
    parser.add_argument("--input_cross_ref", type=str, help="The cross reference table containing contig id, taxonomic annotation and marker gene information for each contig.")
    parser.add_argument("--input_original_bins", type=str, help="The directory containing the original bins from binny.")
    parser.add_argument("--input_MMseqs2_prediction", type=str, help="The MMseqs2 prediction file for the contigs.")
    parser.add_argument("--output_dir_for_refined_bins", type=str, help="The output directory for the refined bins.")
    parser.add_argument("--output_summary_info", type=str, help="The output file for the summary information of the refinement process.")
    parser.add_argument("--log_file", type=str, help="The log file for recording the refinement process.")
    args = parser.parse_args()

    original_bins = args.input_original_bins + "/bins"
    MMseqs2_prediction = args.input_MMseqs2_prediction + "/out.mmseqs2.tsv"
    Path(args.output_summary_info).touch()
    # Inputs and params.
    logging.basicConfig(
        filename=args.log_file,  # Specify the log file
        filemode='a',                # 'a' for append mode, 'w' for write (overwrite)
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Log format
        level=logging.INFO           # Set the minimum log level to INFO
    )

    logger = logging.getLogger(__name__)

    os.mkdir(args.output_dir_for_refined_bins)

    if (os.stat(args.input_cross_ref).st_size == 0):
        print("Empty cross reference table, exit.\n")
        sys.exit(0)

    cross_ref_refine_df = pd.read_csv(args.input_cross_ref, sep="\t", index_col=0)
    bin_list = listdir(original_bins)
    columns = ["contig id", "db_match", "taxonomic_rank", "name", "num_fragments", "num_frag_assigned", "num_agreeing_fragments", "score", "lineage"]
    orig_MMseqs2_df = pd.read_csv(MMseqs2_prediction, sep="\t", names=columns) # we will assume the run used --tax-lineage 2, so the lineage is in the column "lineage".
    lineages = orig_MMseqs2_df['lineage']

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
            elem_row = orig_MMseqs2_df.loc[orig_MMseqs2_df["contig id"] == elem.id]
            elem_lineage = elem_row['lineage'].values[0] if len(elem_row) > 0 else np.nan
            bacteria_taxonomy_id = BACTERIA_ID
            if not pd.isna(elem_lineage):
                elem_lineage_hierarchy = elem_lineage.split(";")
                logging.info("CAT/BAT taxonomic assignment of contig {}: {}".format(elem_id, elem_lineage_hierarchy))
                logging.info("Marker genes found on the contig: {}".format(np.array(elem_in_crossref['markers on the contig'])[0]))
                if (bacteria_taxonomy_id in elem_lineage_hierarchy) and (pd.isna(np.array(elem_in_crossref['markers on the contig'])[0])):
                    logging.info("contig id: {} This is a bacterial contamination.".format(elem_id))
                    contig_ids.remove(elem_id)
                    to_be_modified_id.append(elem_id)
                    with open(args.output_summary_info, "a") as si:
                        si.write(f"MAG {i.split('/')[-1]} has suspicious contamination contig: {elem_id}. Taxonomy: {elem_lineage_hierarchy}.\n")

                elif (EUKARYA_ID in elem_lineage_hierarchy) and (pd.isna(np.array(elem_in_crossref['markers on the contig'])[0])):
                    logging.info("Contig id: {}. This is a eukaryotic contig without marker genes on it.".format(elem_id))
                    contig_ids.remove(elem_id)
                    to_be_modified_id.append(elem_id)
                    with open(args.output_summary_info, "a") as si:
                        si.write(f"MAG {i.split('/')[-1]} has a eukaryotic contig without single-copy marker genes in database: {elem_id}. Taxonomy: {elem_lineage_hierarchy}.\n")

                elif (ROOT_ID in elem_lineage_hierarchy) and (pd.isna(np.array(elem_in_crossref['markers on the contig'])[0])):
                    logging.info("Contig id: {}. This is a contig with ambiguous taxonomic classification and without recognizable A2K marker genes in selected lineage, discarded.\n".format(elem_id))
                    contig_ids.remove(elem_id)
                    to_be_modified_id.append(elem_id)
                    with open(args.output_summary_info, "a") as si:
                        si.write(f"MAG {i.split('/')[-1]} has an unclassified contig: {elem_id}.\n")

                else:
                    continue
            else:
                if elem_in_crossref["markers on the contig"].isna().item():
                    logging.info("Contig id: {}. This is a contig with no taxonomic annotation and no marker genes on it. We will keep it for now, but it is recommended to check the contig manually.\n".format(elem_id))
                    contig_ids.remove(elem_id)
                    to_be_modified_id.append(elem_id)
                    with open(args.output_summary_info, "a") as si:
                        si.write(f"MAG {i.split('/')[-1]} has a contig with no taxonomic annotation and no marker genes: {elem_id}. This contig is removed from the refined bin.\n")
                else:
                    with open(args.output_summary_info, "a") as si:
                        si.write(f"MAG {i.split('/')[-1]} has a contig with no taxonomic annotation: {elem_id}, but has marker genes.\n")
                    logging.info("Contig id: {}. This is a contig with no taxonomic annotation. We will keep it for now, but it is recommended to check the contig manually.\n".format(elem_id))
        
        logging.info(f"contig composition after refinement: {','.join(contig_ids)}\n")
        logging.info(f"length equals to {len(contig_ids)}\n")
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

        
    for i in refined_contigs_dict.items():
        filename = args.output_dir_for_refined_bins + "/" + "refined_" + i[0].split("/")[-1]    
        with open(filename, "w") as refined:
            for j in i[-1]:
                logging.info(j.id)
                refined.write(">")
                refined.write(j.id)
                refined.write("\n")
                refined.write(str(j.seq))
                refined.write("\n")


    cross_ref_refine_df.to_csv(args.input_cross_ref, sep="\t")
    return

if __name__ == "__main__":
    main()