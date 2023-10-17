#!/bin/bash

# exit script if any command results in an error
set -e

# Get the batch ID from the JSON file from crunch
export BATCH_ID=$(python -c "import json;f = open('.crunch/dataset.json') ; print(json.load(f)['name']) ; f.close() ")

if [ -z "$BATCH_ID" ]; then
    echo "Error reading batch id."
    exit 1
fi

mkdir -p $BATCH_ID


# Download the data from mediaflux
unimelb-mf-download --mf.config ~/.Arcitecta/mflux.cfg --csum-check --nb-workers 4 --out $BATCH_ID \
    /projects/proj-6300_metagenomics_processed_verbruggenmdap-1128.4.295/metaGenPipe_outputs/${BATCH_ID}_outputs

# Run Chloroscan, needs to run MMA. Can be solved later. 
bash /data/MMA_organelle_metagenomics/ChloroScan.sh -a "$BATCH_ID/${BATCH_ID}_outputs/assembly/merged.megahit.contigs.fa" -b "chloroscan_outputs/${BATCH_ID}_outputs" -m "$BATCH_ID/${BATCH_ID}_outputs/readalignment/*.bam" -n 500 -k 1 -p 0.90 -t 12
# snakemake -c 8 --snakefile /data/MMA_organelle_metagenomics/Snakefile --directory $BATCH_ID

# Upload the results to mediaflux
unimelb-mf-upload --mf.config ~/.Arcitecta/mflux.cfg --csum-check --nb-workers 4 --create-namespaces \
    --namespace  /projects/proj-6300_metagenomics_processed_verbruggenmdap-1128.4.295/chloroscan_outputs \
    chloroscan_outputs/${BATCH_ID}_outputs

rm -rf $BATCH_ID
rm -rf chloroscan_outputs/${BATCH_ID}_outputs
