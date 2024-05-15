#! /usr/bin/bash 
# Make sure that one error ends all processes.
set -e

Help()
{
    echo "AutoInitiation syntax: bash autoInit.sh [-h|a|b|c|d|s]"
    echo
    echo "The absolute paths below are recommended to be saved first in a text file for convenience."
    echo "-a     Absolute path to your ChloroScan directory."
    echo "-b     Absolute path to your binny checkm_pf hmm profile file."
    echo "-c     Absolute path to your taxon marker set tsv file."
    echo "-d     Database's path where your CAT database is installed or where to install it."
    echo "-h     Print this help message."
    echo "-s     Snakemake virtual environment if you want to specify, for binny installation."
}

# Make sure custom paths to any directories should be added below in terms of ABSOLUTE PATH.
while getopts ":ha:b:c:d:s:" option
do
    case $option in
        h) 
            Help
            exit;;
        a)
            ABS_CHLOROSCAN_DIR=$OPTARG;;
        b)
            PATH_TO_YOUR_BINNY_HMM_PF_DB=$OPTARG;;
        c)
            PATH_TO_YOUR_TAXON_MARKER_SET_TSV=$OPTARG;;
        d)
            DB_DIR=$OPTARG;;
        s)
            snakemake_env=$OPTARG;;
    esac
done

echo "$ABS_CHLOROSCAN_DIR"
echo "$PATH_TO_YOUR_BINNY_HMM_PF_DB"
echo "$PATH_TO_YOUR_TAXON_MARKER_SET_TSV"
echo "$DB_DIR"
echo $snakemake_env

# 1. you have already cloned the raw ChloroScan workflow from the github page, now you should set up all tools.
CHLOROSCAN_DIR="/home/student.unimelb.edu.au/yuhtong/andy/MMA_organelle_metagenomics"
# CHLOROSCAN_DIR="/home/student.unimelb.edu.au/yuhtong/andy/chloroscan"

echo "Current working directory"
pwd

# So now, clone binny into the place.
# make it an if statement.
BINNY_DIR="binny_Chloroscan"
if [ -d $BINNY_DIR ]; then
    echo "Binny for ChloroScan already exists, no need to transfer from github again."
else
    git clone https://github.com/Andyargueasae/binny_Chloroscan.git
fi

cd binny_Chloroscan

# binny initialization should be done.
sed -i "s@snakemake_env: \"\"@snakemake_env: \"$snakemake_env\"@g" config/config.init.yaml

# Check if the conda envs are all there.
if [ ! -z "$(ls -A ./conda)" ]
then
    echo "Binny for ChloroScan already set up, no need to reinstall packages again."
else 
    ./binny -i config/config.init.yaml
fi

# Now move to the change of binny_mantis.cfg and its database.
# first, make sure that mantis only uses checkm_pf as the only hmm database.
echo
echo "Now change the config setup for marker gene usage of binny."
CURR_WORKING_DIR=$(pwd)

echo $CURR_WORKING_DIR
TF_DB_DIR="./database/hmms/checkm_tf/checkm_filtered_tf.hmm"
CFG_PATH="./config/binny_mantis.cfg"
SAMPLE_CFG="sample.cfg"
sed -i "s|custom_ref=$CURR_WORKING_DIR/$TF_DB_DIR||g" $CFG_PATH
sed -i "s|checkm_filtered_tf_weight=0.5||g" $CFG_PATH
sed -i '/^$/d' $CFG_PATH 
# config sorted, now should change the database. Change the hmm profile for pf only.

# You should ensure there is hmmtools in your envs, and only move checkm_pf.hmm in.
BINNY_DB_PATH="./database/hmms/checkm_pf"
# Should make this hmm to the github page.

echo "Are we here?"
rm $BINNY_DB_PATH/checkm_filtered_pf.hmm
# HMM_PROFILE="/data/binny_last_db/hmms/checkm_pf/checkm_filtered_pf.hmm" 
cp $PATH_TO_YOUR_BINNY_HMM_PF_DB $BINNY_DB_PATH
cd $BINNY_DB_PATH/chunks && rm *
cd ../
hmmpress -f checkm_filtered_pf.hmm
cp checkm_filtered_pf.hmm chunks/checkm_filtered_pf_chunk_0.hmm && hmmpress -f chunks/checkm_filtered_pf_chunk_0.hmm

cd ../../
cp $PATH_TO_YOUR_TAXON_MARKER_SET_TSV .

cd $ABS_CHLOROSCAN_DIR
# 2. Install curl and hence frag_gene_scan_rs can be loaded.
#First let's initialize cargo and FragGeneScanRs, fix it to PATH.
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# Should export to PATH.

. "$HOME/.cargo/env"
cargo install frag_gene_scan_rs

CHLOROSCAN_DB=$ABS_CHLOROSCAN_DIR/databases

# Get started by having your data installed.
if [ -n "ls -A $CHLOROSCAN_DB" ]
then
    echo "Your database is already setup."
else
    DB_DIR="/data"
    cd $DB_DIR
    mkdir -p $CHLOROSCAN_DB
    if [ -d "20231120_CAT_nr" ]
    then
        echo "The database has already been installed, you can specify its path and directly use it."
        cd $CHLOROSCAN_DB
        ln -s $DB_DIR/20231120_CAT_nr 20231120_CAT_nrDB

    else
        wget tbb.bio.uu.nl/tina/CAT_pack_prepare/20231120_CAT_nr.tar.gz
        tar -xvzf 20231120_CAT_nr.tar.gz
        rm 20231120_CAT_nr.tar.gz
        # save a bit place.
        cd $CHLOROSCAN_DB
        ln -s $DB_DIR/20231120_CAT_nr 20231120_CAT_nrDB
    fi

fi

cd $ABS_CHLOROSCAN_DIR

# Before you go, check if there are envs called "kronatools".
if [ -n "mamba env list | grep kronatools" ]
then
    echo "There is an virtual environment named kronatools."
else
    echo "Automatically creating the env for kronatools."
    mamba create -n kronatools -c bioconda -c conda-forge -y krona=2.8.1 numpy pandas
fi

# Finally, set up ChloroScan workflow envs.
if [ -n "ls -A ./.snakemake/conda" ]
then
    echo "ChloroScan workflow's envs are already set up, you can directly use it."
else
    echo "Initiate the workflow virtual environments."
    snakemake -c 11 --use-conda --conda-prefix ./conda --snakefile ./chloroscan/workflow/Snakefile --configfile config/ChloroScan.init.yaml --conda-create-envs-only --dryrun
    echo "Workflow envs initiation done. "
fi
# Phew, we have done working out building all envs. Now the ChloroScan workflow is ready to be executed. 