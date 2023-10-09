
default_binny_config_path="binny/config/config.default.yaml"

pwd
# make sure there are no empty spaces when assigning variables.
echo "The default config file path: $default_binny_config_path"

usage() {
    echo "Usage: $0 -[assem|align|output|hdb_epsilon|hdb_min_sample|min_comp|start_comp|purity] " 1>&2
    echo "       -a specify assembly from the config file to be put into the binny config" 1>&2
    echo "       -l specify alignment bam files" 1>&2
    echo "       -o output directory for binny to output those bins" 1>&2
    echo "       -e specify an epsilon range for binny arguments" 1>&2
    echo "       -m specify minimum sample list for each round" 1>&2
    echo "       -c minimum completeness" 1>&2
    echo "       -s start completeness" 1>&2
    echo "       -p purity threshold for bin selection" 1>&2
    echo "       -u universival cutoff" 1>&2
}

while getopts a:d:l:o:e:m:c:s:p:u: flag
do
    case "${flag}" in
        a) assembly=${OPTARG};;
        d) default_config=${OPTARG};;
        l) alignment=${OPTARG};;
        o) outputdir=${OPTARG};;
        e) hdb_epsilon=${OPTARG};;
        m) hdb_min_sample=${OPTARG};;
        c) min_completeness=${OPTARG};;
        s) start_completeness=${OPTARG};;
        p) purity=${OPTARG};;
        u) lengthcutoff=${OPTARG};;
    esac
done

echo "assembly points to: $assembly";
echo "offered alignment files: $alignment";
echo "output directory: $outputdir";
echo "hdbscan epsilon range: $hdb_epsilon";
echo "hdbscan min sample: $hdb_min_sample";
echo "minimum completeness of bin: $min_completeness";
echo "start completeness of bin: $start_completeness";
echo "purity chosen: $purity";
echo "Universal Length cutoff: $lengthcutoff"

# Now edit the file one by one.
## copy that default file and change the copy to be the one.
curr_config="config.ChloroScan.yaml"
# error comes from here first.
echo
echo
echo -------------------------------------------
pwd
echo $default_config
cp $default_config $curr_config
#curr_config=$default_config
# change assembly.
# echo "assembly: \"\""
# create a variable representing "".
double_quotation=\"\"

sed -i "s@assembly: \"\"@assembly: \"$assembly\"@g" $curr_config

#Note this time, binny is nested within one workflow, so it is one folder deeper.
sed -i "s@metagenomics_alignment: \"\"@metagenomics_alignment: \"$alignment\"@g" $curr_config

# Note: the binny is nested within the MMA, if you want to save it inside MMA's output dir, you shall skip back one folder.
sed -i "s@outputdir: \"output\"@outputdir: \"$outputdir\"@g" $curr_config

sed -i "s@min_cont_length_cutoff: 2250@min_cont_length_cutoff: $lengthcutoff@g" $curr_config
awk '/min_cont_length_cutoff/' $curr_config

sed -i "s@max_cont_length_cutoff: 2250@max_cont_length_cutoff: $lengthcutoff@g" $curr_config
awk '/max_cont_length_cutoff/' $curr_config

sed -i "s@min_cont_length_cutoff_marker: 2250@min_cont_length_cutoff_marker: $lengthcutoff@g" $curr_config
awk '/min_cont_length_cutoff_marker/' $curr_config

sed -i "s@max_cont_length_cutoff_marker: 2250@max_cont_length_cutoff_marker: $lengthcutoff@g" $curr_config
awk '/max_cont_length_cutoff_marker/' $curr_config

sed -i "s@hdbscan_epsilon_range: '0.250,0.000'@hdbscan_epsilon_range: '$hdb_epsilon'@g" $curr_config

sed -i "s@hdbscan_min_samples_range: '1,5,10'@hdbscan_min_samples_range: '$hdb_min_sample'@g" $curr_config

sed -i "s@min_completeness: 72.5@min_completeness: $min_completeness@g" $curr_config

sed -i "s@start_completeness: 92.5@start_completeness: $start_completeness@g" $curr_config

sed -i "s@purity: 95@purity: $purity@g" $curr_config

sed -i "s@include_depth_initial: 'False'@include_depth_initial: 'True'@g" $curr_config

sed -i "s@include_depth_main: 'False'@include_depth_main: 'True'@g" $curr_config
# by applying \ before any characters, could help to identify them with efficiency.

# Work done. Now the yaml file is ready to be put into binny. 
# Note: must specify everything in rule in terms of params and input.

# We have to move it back to binny's config folder in order to let it run.
#mv $curr_config $UNIVERSAL_CONFIG_FOR_MMA

# Now you should see it in your binny workflow config directory.