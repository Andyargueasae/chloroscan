#!env bash
set -o pipefail
set -e 

git config --global http.postBuffer 524288000

git clone https://github.com/Andyargueasae/binny_Chloroscan.git $CONDA_PREFIX/lib/binny_Chloroscan

CHLOROSCAN_DIR=$(pwd)

cd $CONDA_PREFIX/lib/binny_Chloroscan

BINNY_DIR=$(pwd)

./binny -i config/config.init.yaml

TF_DB_DIR="database/hmms/checkm_tf/checkm_filtered_tf.hmm"
CFG_PATH="./config/binny_mantis.cfg"
SAMPLE_CFG="sample.cfg"
sed -i "s|custom_ref=$BINNY_DIR/$TF_DB_DIR||g" $CFG_PATH
sed -i "s|checkm_filtered_tf_weight=0.5||g" $CFG_PATH
sed -i '/^$/d' $CFG_PATH 
# config sorted, now should change the database. Change the hmm profile for pf only.

REPLACEMENT_HMM=$CHLOROSCAN_DIR/binny_algal_database/checkm_filtered_pf.hmm

CHECKM_PF_DIR="database/hmms/checkm_pf"
rm ./database/hmms/checkm_pf/checkm_filtered_pf.hmm 
cp $REPLACEMENT_HMM $CHECKM_PF_DIR
hmmpress -f $CHECKM_PF_DIR/checkm_filtered_pf.hmm

rm -rf $CHECKM_PF_DIR/chunks/*
cp $CHECKM_PF_DIR/checkm_filtered_pf.hmm $CHECKM_PF_DIR/chunks/checkm_filtered_pf_chunk_0.hmm && hmmpress -f $CHECKM_PF_DIR/chunks/checkm_filtered_pf_chunk_0.hmm

REPLACEMENT_TSV=$CHLOROSCAN_DIR/binny_algal_database/taxon_marker_sets_lineage_sorted.tsv
cp $REPLACEMENT_TSV ./database

# export PATH=$BINNY_DIR:$PATH
 
# activate script PATH
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
echo -e "export BINNY_DIR=\"$BINNY_DIR\"" > "$CONDA_PREFIX/etc/conda/activate.d/binny.sh"


