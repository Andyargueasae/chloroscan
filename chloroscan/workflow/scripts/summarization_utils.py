from Bio import SeqIO

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
        contig, depth = depth_lines[i].strip().split("\t",1)[0], depth_lines[i].strip().split("\t",1)[1]
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
            the_organism = all_organisms[id_hierarchy[i]]
            i+=1
        except KeyError:
            break
    return all_organisms[id_hierarchy[i-1]]

def all_organisms_dict(names_dump):
    all_organisms = dict()
    with open(names_dump, "r") as nd:
        org_lines = nd.readlines()
    
    for i in org_lines:
        i = i.strip().replace("\t", "")
        org_line = i.split("|")
        if "scientific name" in org_line:
            all_organisms[org_line[0]] = org_line[1]
    return all_organisms

def form_contig_taxa_array(empty_taxon_array, contig_id, names_dump, contig_annotation):
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