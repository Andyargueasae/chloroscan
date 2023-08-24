# -*- coding: utf-8 -*-
"""
plastid contigs filter.
"""

import Bio
from Bio import SeqIO
import numpy as np
import os
from os import listdir
import sys
import pandas as pd
import click

@click.command()
#Get all parameters. Especially the depth file.

@click.argument('contigs', type = click.Path(exists=True))
@click.argument('corgipred', type = click.Path(exists=True))
@click.option("--pthreshold", type = click.FloatRange(0, 1), default = 0.25)
@click.option("--filenumb", default = 2, type = int)
@click.option("--output", default = "plastid_contigs_mk2.fasta", type=str)
@click.option("--depthfile", default = None, type = click.Path(exists = True))
@click.option("--depthout", default = "new_depth.txt", type=str)

def main(contigs, corgipred, pthreshold, filenumb, output, depthfile, depthout):
    """
    Collect the contigs predicted to be plastidial and store in a new fasta 
    file.
    """
    click.echo(f"contigs file: {contigs}")
    click.echo(f"corgi prediction passed as: {corgipred}")
    click.echo(f"There are {filenumb} files specified.\n")
    click.echo(f"The threshold for a contig to be plastidial is: {pthreshold}")
    
    click.echo("Now initiate the plastid contigs filtration.")
    
    #Feed the SeqIO with final contigs specified in the part.
    final_contigs = list(SeqIO.parse(handle=contigs, format = "fasta"))
    #The csv prediction output by corgi is also loaded, now we're good to go.
    corgi_prediction = pd.read_csv(corgipred, header = 0)
    plastid_id = plastid_filter(final_contigs, corgi_prediction, pthreshold, 
                                output)
    print(f"The process is done. File name: {output}")
    if (depthfile != None) and (depthout):
        click.echo(f"depth file generated from metabat2:{depthfile}")
        depth_filter(depthfile, plastid_id, depthout)
        click.echo(f"depth file updated in {depthout}!")
    return

def plastid_filter(contigs, predictions, p_thresh, output):
    # This function does all the things: filter out plastid contigs and store
    # them in a fasta file.
    contig_id=[]
    contig_seq=[]
    for i in contigs:
        contig_id.append(i.id)
        contig_seq.append(i.seq)
    
    #Now, filter those are plastidial.
    candidate_plastidial = predictions.loc[predictions['prediction'] == "plastid"]
    #print(candidate_plastidial)
    candidate_plastid_contigs = list(candidate_plastidial['accession'])
    #print(candidate_plastid_contigs)
    length_cand = len(candidate_plastid_contigs)
    print(f"Found {length_cand} candidate plastidial contigs.")
    plastid_contig_seqs = []
    count = 0
    plastid_contigs_id = []
    for elem in candidate_plastid_contigs:
        index = contig_id.index(elem)
        plastid_contigs_id.append(index)
        plastid_contig_seqs.append(contig_seq[index])
        
        count += 1
    
    #Now at here we've done storing sequences along with id. Write a file.
    plastid_fasta_mk2 = open(output, "w")
    for i in range(len(candidate_plastid_contigs)):
        plastid_fasta_mk2.write(">" + candidate_plastid_contigs[i] +
                                "\n" + str(plastid_contig_seqs[i]) + "\n")
    plastid_fasta_mk2.close()
    # May return contig_id as part of response.
    return plastid_contigs_id

def depth_filter(depthfile, plastid_id_list, depthout):
    '''Get the depth.txt file as soon as possible, maybe tomorrow!'''
    print("Passing in contig id in order of plastid contigs.")
    id_array = np.array(plastid_id_list)
    with open(depthfile) as dpf:
        lines_depth = dpf.readlines()
    #Headers shall be restored.
    header = lines_depth[0]
    lines_depth.remove(header)
    
    print("Now the lines are located and written into the file.")
    filtered_depth = open(depthout, "w")
    filtered_depth.write(header)
    for i in id_array:
        line_to_write = lines_depth[i]
        filtered_depth.write(line_to_write)
    filtered_depth.close()
        
    #New depth_file is ready.
    return
    
if __name__ == "__main__":
    main()
