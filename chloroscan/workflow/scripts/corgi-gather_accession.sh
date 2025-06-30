#!/bin/bash

# This bash script takes in two steps:
# 1. takes in a csv file as the input to extract those contig id with a prediction to be plastidial.
# 2. take those contig ids to write a fasta file, containing all plastidial contigs.

# This script can be transferred to either the inputing-style or snakemake script style.

plastid_key="plastid"
bacteria_control=0.10

while getopts c:o:p: flag
do
    case "${flag}" in
        c) prediction_csv=${OPTARG};;
        o) output_path=${OPTARG};;
        p) pthreshold=${OPTARG};;
    esac
done

while IFS=, read -r file accession archaea bacteria eukaryote mitochondrion plastid prediction eukaryotic prokaryotic organellar
do 
    if [ "$prediction" = "$plastid_key" ] && (( $(echo $plastid $pthreshold | awk '{if ($1 >= $2) print 1;}') )); then
        echo "$accession"
    fi
done < <(tail -n +2 $prediction_csv) > $output_path

