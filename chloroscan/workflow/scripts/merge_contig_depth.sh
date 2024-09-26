# Merge the contig depth.
start=`date +%s`
set -e
first_file=true
input=${snakemake_input}
output=${snakemake_output}
for file in $(ls $input); do
    if [[ $first_file == 'true' ]]; then
        # echo "First file."
        cp $input/$file $output
        first_file=false
    else
        # echo "File $COUNTER"
        export TMPDIR="MergeDepthTMPDIR"
        mkdir -p $TMPDIR
        TMP_DEPTH=$(mktemp --tmpdir=$TMPDIR "tmp_XXXXXXXXXXX.tsv")
        paste $output <( cut -f 2 $input/$file) > $TMP_DEPTH \
        && mv $TMP_DEPTH $output
    fi
done

end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo $runtime
