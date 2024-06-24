#! usr/bin/bash
set -e
default_binny_config_path="binny_ChloroScan/config/config.default.yaml"

pwd
# make sure there are no empty spaces when assigning variables.
echo "The default config file path: $default_binny_config_path"

usage() {
    echo "Usage: $0 -[h|a|d|t|l|o|e|m|c|s|p|u|y] " 1>&2
    echo "       -h Show this help message." 1>&2
    echo "       -a specify assembly from the config file to be put into the binny config" 1>&2
    echo "       -b specify pre-downloaded binny database to use." 1>&2
    echo "       -d (Absolute) Path to default binny configuration used for template" 1>&2
    echo "       -t Depth profile for assembly" 1>&2
    echo "       -l specify alignment bam files" 1>&2
    echo "       -o output directory for binny to output those bins" 1>&2
    echo "       -e specify an epsilon range for binny arguments" 1>&2
    echo "       -m specify minimum sample list for each round" 1>&2
    echo "       -c minimum completeness" 1>&2
    echo "       -s start completeness" 1>&2
    echo "       -p purity threshold for bin selection" 1>&2
    echo "       -u universival cutoff" 1>&2
    echo "       -y (Absolute) Path (and names) of your binny config to be executed." 1>&2
}

while getopts ":ha:b:d:t:l:o:e:m:c:s:p:u:y:" option
do
    case $option in
        h) 
            usage
            exit;; 
        a) 
            assembly=$OPTARG;;
        b)
            binny_db=$OPTARG;;
        d) 
            default_config=$OPTARG;;
        t) 
            depth_text=$OPTARG;;
        l) 
            alignment=$OPTARG;;
        o) 
            outputdir=$OPTARG;;
        e) 
            hdb_epsilon=$OPTARG;;
        m) 
            hdb_min_sample=$OPTARG;;
        c)
            min_completeness=$OPTARG;;
        s) 
            start_completeness=$OPTARG;;
        p) 
            purity=$OPTARG;;
        u) 
            lengthcutoff=$OPTARG;;
        y)
            curr_config=$OPTARG;;
    esac
done

echo "assembly points to: $assembly";
echo "Default binny config file: $default_config";
echo "Config file to be used: $curr_config";
echo "offered alignment files: $alignment";
echo "offered depth profile: $depth_text"
echo "output directory: $outputdir";
echo "hdbscan epsilon range: $hdb_epsilon";
echo "hdbscan min sample: $hdb_min_sample";
echo "minimum completeness of bin: $min_completeness";
echo "start completeness of bin: $start_completeness";
echo "purity chosen: $purity";
echo "Universal Length cutoff: $lengthcutoff"

# Now edit the file one by one.
## copy that default file and change the copy to be the one.
# error comes from here first.
echo
echo -------------------------------------------
pwd
cp $default_config $curr_config
double_quotation=\"\"

# Must make sure the config is going to back to its directory. 

echo "$depth_text" # The files are expanded.
echo "$alignment"
# The script assumes that either alignment bams or depth text files are supplied.
if [ -n "$depth_text" ] && [ -z "$alignment" ];
then
    echo "The alignment bam file(s) have already been transferred into tab-separated text file, use txt."
    sed -i "s@contig_depth: \"\"@contig_depth: \"$depth_text\"@g" $curr_config
elif [ -n "$alignment" ] && [ -z "$depth_text" ];
then
    echo "Only alignment bam file(s) have been provided."
    sed -i "s@metagenomics_alignment: \"\"@metagenomics_alignment: \"$alignment\"@g" $curr_config
elif [ -n "$alignment" ] && [ -n "$depth_text" ];
then 
    echo "Both alignment bam files and depth profile txt are provided, use depth profile for downstream jobs."
    sed -i "s@contig_depth: \"\"@contig_depth: \"$depth_file\"@g" $curr_config
else
    echo "No alignment bam files or depth profiles were provided, binny cannot proceed."
    exit 0
fi

# Change the db path.
sed -i "s@db_path: \"\"@db_path: \"$binny_db\"@g" $curr_config

sed -i "s@assembly: \"\"@assembly: \"$assembly\"@g" $curr_config

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

# Now you should see it in your binny workflow config directory.