#!/usr/bin/env bash

# The binny directory.
BINNY_DIR=${snakemake_input[binny_dir]}/bins

# assign variables.
FILE_FLAG=${snakemake_params[gff_file_flag]}
BATCH_NAME=${snakemake_params[batch_name]}
# PATH_NAME_BINS=${snakemake_params[path_name_bins]}
FGSR_OUTPUT=${snakemake_params[path_fraggenescanrs]}

CDS_OUTPUT=${snakemake_output[CDS_EXTRACTION]}

mkdir -p $FGSR_OUTPUT/$FILE_FLAG
mkdir -p $CDS_OUTPUT/cds
mkdir -p $CDS_OUTPUT/faa
echo ${snakemake_output}


for i in $(ls $BINNY_DIR);
do 
    echo $BINNY_DIR/$i
    # echo ${snakemake_output}
    FragGeneScanRs --seq-file-name $BINNY_DIR/$i --training-file illumina_5 --thread-num 1 -g $FGSR_OUTPUT/$FILE_FLAG/$(basename $i .fasta).gff \
     -n $CDS_OUTPUT/cds/$(basename $BATCH_NAME).$(basename $i .fasta).gene.fasta -a $CDS_OUTPUT/faa/$(basename $BATCH_NAME).$(basename $i .fasta).faa -w 0; 
done
