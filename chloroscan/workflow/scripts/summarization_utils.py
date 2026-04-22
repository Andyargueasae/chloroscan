from Bio import SeqIO
import pandas as pd
import numpy as np

def GC_content(sequence):
    G_count = sequence.count('G')
    C_count = sequence.count('C')
    return (G_count + C_count)/len(sequence)

def form_bin_array(bin_dir, collected_bins, bin_array, contig_id_list):
    """
    bin_dir: bin_directory inside binny;
    collected_bins: bin_file name;
    bin_array: the empty list containing len(contig_id) elements.
    contig_id_list: the list of contigs with their id.
    """
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
    with open(annot_file) as af:
        annot_lines = af.readlines()
    marker_dict = {}
    default = "no markers in the database"
    for i in range(len(annot_lines)):
        i_th_ORF = annot_lines[i].split("\t")
        
        key = i_th_ORF[0]
        ORF = i_th_ORF[-1]
        ORF_split = ORF.strip().split(";")
        # Initiate the gathering.
        if ("checkm_marker" in ORF_split[-1]) and (key not in marker_dict.keys()):
            marker_dict[key] = ""
            marker_dict[key] += ORF_split[-1].split("=")[-1]
        elif ("checkm_marker" in ORF_split[-1]) and (key in marker_dict.keys()):
            marker_dict[key] += ORF_split[-1].split("=")[-1]
            marker_dict[key] += ","
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

    return empty_marker_array

def form_depth_array(empty_depth_array, contig_id, depth_file):
    with open(depth_file) as df:
        depth_lines = df.readlines()
    for i in range(len(depth_lines)):
        contig, depth = depth_lines[i].strip().split("\t",1)[0], depth_lines[i].strip().split("\t",1)[1]
        if contig in contig_id:
            empty_depth_array[contig_id.index(contig)] = depth
    return empty_depth_array

def store_basic_info(contig_id_array, contig_seq_array, contig_len_array, contig_gc_array, assembly_contigs):
    """
    contig_id_array: the empty list containing len(contig_id) elements. To save the contig id.
    contig_seq_array: the empty list containing len(contig_id) elements. To save the contig sequence.
    contig_len_array: the empty list containing len(contig_id) elements. To save the contig length.
    contig_gc_array: the empty list containing len(contig_id) elements. To save the contig GC content.
    assembly_contigs: the contig dictionary.
    """
    # simply change all of the files
    for i in assembly_contigs:
        contig_id_array.append(i.id)
        contig_seq_array.append(i.seq)
        contig_len_array.append(len(i.seq))
        contig_gc_array.append(GC_content(str(i.seq)))
    return

def form_contig_taxa_array(empty_taxon_array: list, contig_id: list, contig_annotation: str):
    """
    empty_taxon_array: the empty list containing len(contig_id) elements.
    contig_id: the list of contig id.
    contig_annotation: the file containing contig id and their corresponding taxid.
    """
    columns = ["contig id", "db_match", "taxonomic_rank", "name", "num_fragments", "num_frag_assigned", "num_agreeing_fragments", "score", "lineage"]
    mmseqs2_contig_taxonomy = pd.read_csv(contig_annotation, sep="\t", names=columns)
    # fill in empty_taxon_array in the order of contig_id list.
    for i in contig_id:
        taxon_row = mmseqs2_contig_taxonomy.loc[mmseqs2_contig_taxonomy["contig id"] == i]
        taxon_id = np.array(taxon_row["name"])[0]
        empty_taxon_array[contig_id.index(i)] = taxon_id
    return