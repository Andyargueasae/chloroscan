#!/bin/bash

# This bash script takes in two steps:
# 1. takes in a csv file as the input to extract those contig id with a prediction to be plastidial.
# 2. take those contig ids to write a fasta file  containing all plastidial contigs.

# This script can be transferred to either the inputing-style or snakemake script style.

# #prediction="plastid"
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--csv)
            prediction_csv="$2"
            shift # past argument
            shift # past value
            ;;
        -o|--output)
            output_path="$2"
            shift # past argument
            shift # past value
            ;;
        -p|--pthreshold)
            pthreshold="$2"
            shift # past argument
            shift # past value
            ;;
        *) # unknown option
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

#file accession prediction probability original_id description Nuclear Mitochondrion Plastid Plasmid Nuclear/Bacteria Nuclear/Archaea Nuclear/Eukaryota Nuclear/Viruses Mitochondrion/Eukaryota Plastid/Eukaryota Plasmid/Bacteria Plasmid/Archaea Plasmid/Eukaryota
plastid_key="Plastid"

# New pthreshold shall be 0.10 to make a trade off between false positive and true positive.

while IFS=, read -r file accession prediction probability original_id description Nuclear Mitochondrion Plastid Plasmid Nuclear_Bacteria Nuclear_Archaea Nuclear_Eukaryota Nuclear_Viruses Mitochondrion_Eukaryota Plastid_Eukaryota Plasmid_Bacteria Plasmid_Archaea Plasmid_Eukaryota
do 
    if [ "$prediction" = "$plastid_key" ] && (( $(echo $Plastid $pthreshold | awk '{if ($1 >= $2) print 1;}') )); then
        echo "$accession"
    fi
done < <(tail -n +2 $prediction_csv) > $output_path
