#!/usr/bin/env bash

# The binny directory.
BINNY_DIR=${snakemake_input[binny_dir]}

# assign variables.
FILE_FLAG=${snakemake_params[gff_file_flag]}
BATCH_NAME=${snakemake_params[batch_name]}
# PATH_NAME_BINS=${snakemake_params[path_name_bins]}
FGSR_OUTPUT=${snakemake_params[path_fraggenescanrs]}

CDS_OUTPUT=${snakemake_output[CDS_EXTRACTION]}
END_FLAG=${snakemake_output[END_FLAG]}
mkdir -p $FGSR_OUTPUT/$FILE_FLAG

# iterate over each bin recovered by binny.
for i in $(ls $BINNY_DIR);
do 
    echo $BINNY_DIR/$i
    # echo ${snakemake_output}
    FragGeneScanRs --seq-file-name $BINNY_DIR/$i --training-file illumina_5 --thread-num 1 -g $FGSR_OUTPUT/$FILE_FLAG/$(basename $i .fasta).gff -w 0; 
done

# rm $FGSR_OUTPUT/$FILE_FLAG/*.fasta.fai.gff
# rm $BINNY_DIR/$PATH_NAME_BINS/*.fai
# Don't forget the end step, change the permission of all gff outputs, as if you don't make the change they will not be able to opened by gffread.
# This script could make use of snakemake parameters to give cds extraction for those MAGs.

# make the directory first.
mkdir -p $CDS_OUTPUT/cds
echo ${snakemake_output}

for i in $(ls $FGSR_OUTPUT/$FILE_FLAG);
do
    echo $i
    # gffread -C ${snakemake_input}/$i -o $(basename $i .gff).cds 2> ${snakemake_log};
    gffread -w $(basename $BATCH_NAME).$(basename $i .gff).cds.fasta -g $BINNY_DIR/$(basename $i .gff).fasta $FGSR_OUTPUT/$FILE_FLAG/$i
    mv $(basename $BATCH_NAME).$(basename $i .gff).cds.fasta $CDS_OUTPUT/cds
    rm $BINNY_DIR/*.fai

done

# So basically, every script uses a for loop to do batch jobs for each of the files.
## As this rule is the last rule of the workflow, touch one final file within output.
# We can't just let the flag file come out, we must need to test whether there are files in cds here!
touch $END_FLAG