#! /usr/bin/bash 

# Make sure custom paths to any directories should be added below in terms of ABSOLUTE PATH.
while getopts: "a:b:c:d" flag
do
    case "${flag}" in
        a) ABS_CHLOROSCAN_DIR=${OPTARG};;
        b) PATH_TO_YOUR_BINNY_HMM_PF_DB=${OPTARG};;
        c) PATH_TO_YOUR_TAXON_MARKER_SET_TSV=${OPTARG};;
        d) DB_DIR=${OPTARG};;
        r) CHLOROSCAN_DIR=${OPTARG};;
    esac
done

echo "$ABS_CHLOROSCAN_DIR"
echo "$PATH_TO_YOUR_BINNY_HMM_PF_DB"
echo "$PATH_TO_YOUR_TAXON_MARKER_SET_TSV"
echo "$DB_DIR"
echo "$CHLOROSCAN_DIR"

# 1. you have already cloned the raw ChloroScan workflow from the github page, now you should set up all tools.
CHLOROSCAN_DIR="/home/student.unimelb.edu.au/yuhtong/andy/MMA_organelle_metagenomics"
# CHLOROSCAN_DIR="/home/student.unimelb.edu.au/yuhtong/andy/chloroscan"

# So now, clone binny into the place.
# make it an if statement.
BINNY_DIR="binny_Chloroscan"
if [-d $BINNY_DIR]; then
    echo "Binny for ChloroScan already exists, no need to transfer from github again."
else
    git clone https://github.com/Andyargueasae/binny_Chloroscan.git
fi

cd binny_Chloroscan

# binny initialization should be done.
./binny -i config/config.init.yaml
# Now move to the change of binny_mantis.cfg and its database.
# first, make sure that mantis only uses checkm_pf as the only hmm database.
CURR_WORKING_DIR=$(pwd)
DB_DIR="binny_ChloroScan/database/hmms/checkm_tf/checkm_filtered_tf.hmm"
CFG_DIR="./config/binny_mantis.cfg"
SAMPLE_CFG="sample.cfg"
sed -i -e "s|custom_ref=$CURR_WORKING_DIR/$DB_DIR||g" $SAMPLE_CFG
sed -i -e "s|checkm_filtered_tf_weight=0.5||g" $SAMPLE_CFG
sed -i '/^$/d' $SAMPLE_CFG 
# config sorted, now should change the database. Change the hmm profile for pf only.

# You should ensure there is hmmtools in your envs, and only move checkm_pf.hmm in.
BINNY_DB_PATH="./binny_ChloroScan/database/hmms/checkm_pf"
Should make this hmm to the github page.
rm $BINNY_DB_PATH/checkm_filtered_pf.hmm
HMM_PROFILE="/data/binny_last_db/hmms/checkm_pf/checkm_filtered_pf.hmm" 
cp $HMM_PROFILE $BINNY_DB_PATH
cd $BINNY_DB_PATH/chunks && rm *
cd ../
hmmpress -f checkm_filtered_pf.hmm
cp checkm_filtered_pf.hmm chunks/checkm_filtered_pf_chunk_0.hmm && hmmpress -f chunks/checkm_filtered_pf_chunk_0.hmm

cd ../../
cp $PATH_TO_YOUR_TAXON_MARKER_SET_TSV .

cd $CHLOROSCAN_DIR
# 2. Install curl and hence frag_gene_scan_rs can be loaded.
#First let's initialize cargo and FragGeneScanRs, fix it to PATH.
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# Should export to PATH.
. "$HOME/.cargo/env"
cargo install frag_gene_scan_rs

# Get started by having your data installed.
DB_DIR="/data"
cd $DB_DIR

CHLOROSCAN_DB=$CHLOROSCAN_DIR/databases
mkdir -p $CHLOROSCAN_DB
if [ -d "20231120_CAT_nr" ]
then
    echo "The database has already been installed, you can specify its path and directly use it."
    ln -s 20231120_CAT_nr $CHLOROSCAN_DB

else
    wget tbb.bio.uu.nl/tina/CAT_pack_prepare/20231120_CAT_nr.tar.gz
    tar -xvzf 20231120_CAT_nr.tar.gz
    rm 20231120_CAT_nr.tar.gz
    # save a bit place.
fi

cd $CHLOROSCAN_DIR
# Finally, set up ChloroScan workflow envs.
snakemake -c 11 --use-conda --snakefile chloroscan/workflow/Snakefile --configfile config/ChloroScan.init.yaml --conda-create-envs-only
# Phew, we have done working out building all envs. Now the ChloroScan workflow is ready to be executed. 