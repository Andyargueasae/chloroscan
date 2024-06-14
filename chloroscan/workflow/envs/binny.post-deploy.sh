#!env bash
set -o pipefail
set -e 

git config --global http.postBuffer 524288000

# We don't even know whether this one is being installed.
git clone https://github.com/Andyargueasae/binny_Chloroscan.git $CONDA_PREFIX/lib/binny_Chloroscan

# cpabc #add a nonsense code to see if cloning was ok.  

CHLOROSCAN_DIR=$(pwd)

cd $CONDA_PREFIX/lib/binny_Chloroscan

BINNY_DIR=$(pwd)

# Now that the code has been changed, we circumvent the download.
# ./binny -i config/config.init.yaml

USE_DB_DIR="binny_algal_database/hmms/checkm_pf/checkm_filtered_pf.hmm"
CFG_PATH="./config/binny_mantis.cfg"
# SAMPLE_CFG="./config/binny_mantis.cfg"
sed -i "s|custom_ref=/mnt/lscratch/users/ohickl/binning/tools/binny_devel/database/hmms/checkm_tf/checkm_filtered_tf.hmm||g" $CFG_PATH
sed -i "s|checkm_filtered_tf_weight=0.5||g" $CFG_PATH
sed -i "s|/mnt/lscratch/users/ohickl/binning/tools/binny_devel/database/hmms/checkm_pf/checkm_filtered_pf.hmm|$BINNY_DIR/$USE_DB_DIR|g" $CFG_PATH
sed -i '/^$/d' $CFG_PATH 
# config sorted, now should setup the database. use mantis only.
mantis setup -mc $CFG_PATH --no_taxonomy 
mantis check -mc $CFG_PATH --no_taxonomy

# configuration file need to be changed. 
# REPLACEMENT_HMM=$BINNY_DIR/binny_algal_database/checkm_filtered_pf.hmm

# CHECKM_PF_DIR="database/hmms/checkm_pf"
# rm ./database/hmms/checkm_pf/checkm_filtered_pf.hmm 
# cp $REPLACEMENT_HMM $CHECKM_PF_DIR
# hmmpress -f $CHECKM_PF_DIR/checkm_filtered_pf.hmm

# rm -rf $CHECKM_PF_DIR/chunks/*
# cp $CHECKM_PF_DIR/checkm_filtered_pf.hmm $CHECKM_PF_DIR/chunks/checkm_filtered_pf_chunk_0.hmm && hmmpress -f $CHECKM_PF_DIR/chunks/checkm_filtered_pf_chunk_0.hmm

# REPLACEMENT_TSV=$BINNY_DIR/binny_algal_database/taxon_marker_sets_lineage_sorted.tsv
# cp $REPLACEMENT_TSV ./database

# export PATH=$BINNY_DIR:$PATH
 
# activate script PATH
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
echo -e "export BINNY_DIR=\"$BINNY_DIR\"" > "$CONDA_PREFIX/etc/conda/activate.d/binny.sh"


