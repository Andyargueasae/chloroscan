#!/usr/bin/env bash
# Now it is time to think about some global variables.
echo ${snakemake_input[bins]}

mkdir -p ${snakemake_output}
for i in $(ls ${snakemake_input[bins]});
do 
    echo ${snakemake_input[bins]}/$i
    echo ${snakemake_output}
    FragGeneScanRs --seq-file-name ${snakemake_input[bins]}/$i --training-file illumina_5 --thread-num 1 -g ${snakemake_output}/$(basename $i .fasta).gff -w 0;
done

