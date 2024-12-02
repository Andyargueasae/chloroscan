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

DEFAULT_BINNY_CONDASOURCE=$BINNY_DIR/conda

default_binny_config_path="$CONDA_PREFIX/lib/binny_Chloroscan/config/config.default.yaml"

# Now that the code has been changed, we circumvent the download.
# ./binny -i config/config.init.yaml
echo "snakemake env: $(which snakemake)"
# Create the envs for snakemake.
mamba create --prefix $BINNY_DIR/snakemake_env -y snakemake=7.16.0 unzip python=3.8 

# Create envs for other binny rules.
# snakemake --use-conda --conda-create-envs-only --conda-prefix $DEFAULT_BINNY_CONDASOURCE --configfile $default_binny_config_path --snakefile $BINNY_DIR/Snakefile --cores 1 --verbose

# Add this new env to config file of binny. 
sed -i "s@snakemake_env: \"\"@snakemake_env: \"$BINNY_DIR/snakemake_env\"@g" $default_binny_config_path


# This is the original code, but we need to change the database to the one A2K one.
USE_DB_DIR="A2K_database/hmms/checkm_pf/checkm_filtered_pf.hmm"
CFG_PATH="./config/binny_mantis.cfg"
# SAMPLE_CFG="./config/binny_mantis.cfg"
sed -i "s|custom_ref=/mnt/lscratch/users/ohickl/binning/tools/binny_devel/database/hmms/checkm_tf/checkm_filtered_tf.hmm||g" $CFG_PATH
sed -i "s|checkm_filtered_tf_weight=0.5||g" $CFG_PATH
sed -i "s|/mnt/lscratch/users/ohickl/binning/tools/binny_devel/database/hmms/checkm_pf/checkm_filtered_pf.hmm|$BINNY_DIR/$USE_DB_DIR|g" $CFG_PATH
sed -i '/^$/d' $CFG_PATH 
# config sorted, now should setup the database. use mantis only.

mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
echo -e "export BINNY_DIR=\"$BINNY_DIR\"" > "$CONDA_PREFIX/etc/conda/activate.d/binny.sh"


