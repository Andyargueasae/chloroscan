#!/usr/bin/env bash
set -euo pipefail
# The binny directory. Is from the parental directory of binny.done, which works as the input. 
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) BINNY_DONE="$2"; shift ;;
        --batch_name) BATCH_NAME="$2"; shift ;;
        --gff_file_flag) FILE_FLAG="$2"; shift ;;
        --path_fraggenescanrs) FGSR_OUTPUT="$2"; shift ;;
        --output) CDS_OUTPUT="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

mkdir -p $FGSR_OUTPUT/$FILE_FLAG
mkdir -p $CDS_OUTPUT/cds
mkdir -p $CDS_OUTPUT/faa
echo ${CDS_OUTPUT}

export PATH=$PATH:$HOME/.cargo/bin

BINNY_DIR=$(dirname $BINNY_DONE)/bins;

for i in $(ls $BINNY_DIR);
do 
    echo $BINNY_DIR/$i
    # echo ${snakemake_output}
    FragGeneScanRs --seq-file-name $BINNY_DIR/$i --training-file illumina_5 --thread-num 1 -g $FGSR_OUTPUT/$FILE_FLAG/$(basename $i .fasta).gff \
     -n $CDS_OUTPUT/cds/$(basename $BATCH_NAME).$(basename $i .fasta).gene.fasta -a $CDS_OUTPUT/faa/$(basename $BATCH_NAME).$(basename $i .fasta).faa -w 0; 
done
