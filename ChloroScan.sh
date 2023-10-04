#!/bin/bash
# This file writes a new configuration yaml file for each batch, and send them into the snakemake workflow.
# Can be customized to write any parameters you want.

set -e

# NOTE: every path should be in relative path.

# A reminder: MUST add quotation for assembly, alignment and batch_name, they may interact with the opts. 
while getopts "a:b:m:n:k:p:t:" flag
do
    case "${flag}" in
        a) assembly=${OPTARG};;
        b) batch_name=${OPTARG};;
        m) alignment=${OPTARG};;
        n) min_length_corgi=${OPTARG};;
        k) corgi_batch_size=${OPTARG};;
        p) p_threshold=${OPTARG};;
        t) threads=${OPTARG};;
    esac
done

# batch_name is used to name the output directory.
mkdir -p config/
cp /data/MMA_organelle_metagenomics/config/config.zero.yaml config/config.default.yaml

input_config="config/config.default.yaml"
default_CAT_DB="CAT_prepare_20210107"
default_taxonomy_dump="names.dmp"


echo "min length is $min_length_corgi"
echo "corgi batch size is $corgi_batch_size"

# Now change some settings in relation to the config file.
# Note: all path should be relative path.
sed -i "s@assembly_path: @assembly_path: \"$assembly\"@g" $input_config
awk "/assembly_path/" $input_config

sed -i "s@alignment_files: @alignment_files: \"$alignment\"@g" $input_config
awk "/alignment_files/" $input_config

sed -i "s@outputdir: \"chloroscan_outputs\"@outputdir: \"$batch_name\"@g" $input_config
awk "/outputdir:/" $input_config

sed -i "s@batch_name: @batch_name: \"$batch_name\"@g" $input_config
awk "/batch_name/" $input_config

sed -i "s@minlen: @minlen: $min_length_corgi@g" $input_config
awk "/minlen/" $input_config

sed -i "s@batch_size: @batch_size: $corgi_batch_size@g" $input_config
awk "/batch_size/" $input_config

# write CAT database to config.

sed -i "s@CAT_database: @CAT_database: \"$default_CAT_DB\"@g" $input_config
awk "/CAT_database/" $input_config

# save taxonomy dump.
sed -i "s@taxonomy_names: @taxonomy_names: \"$default_taxonomy_dump\"@g" $input_config
awk "/taxonomy_names/" $input_config

sed -i "s@pthreshold: @pthreshold: $p_threshold@g" $input_config
awk "/pthreshold/" $input_config

sed -i "s@threads: @threads: $threads@g" $input_config
awk "/threads/" $input_config

# We only need to initiate the snakemake workflow by default.

snakemake -c --latency-wait 15 --use-conda --configfile $input_config \
    --snakefile /data/MMA_organelle_metagenomics/Snakefile \
    --conda-prefix /data/MMA_organelle_metagenomics/conda \
    --conda-frontend mamba

output_dir=$batch_name
cp $input_config $output_dir


# remember always remove it after workflow done. 
rm $input_config