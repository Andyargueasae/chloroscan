start=`date +%s`
set -e
TMPDIR=${snakemake_params[tmpdir]}
# TMP_DEPTH=$(mktemp --tmpdir={TMPDIR} "depth_file_XXXXXXXXXXXXXXXX.txt")

input_assembly=${snakemake_input[assembly]}
OUTPUT_CONTIGDEPTH=${snakemake_output}
[  ! -d $OUTPUT_CONTIGDEPTH ] && mkdir -p $OUTPUT_CONTIGDEPTH
for depth in $(ls ${snakemake_input[alignment]} | grep .bam)
do
    genomeCoverageBed -ibam ${snakemake_input[alignment]}/$depth | grep -v "genome" > $TMPDIR/$(basename $depth .sorted.bam).txt
    echo "Depth calculation got, now perform average contig depth calculation." >> ${snakemake_log}
    perl ${snakemake_params[perl_script]} $TMPDIR/$(basename $depth .sorted.bam).txt $input_assembly > $OUTPUT_CONTIGDEPTH/$(basename $depth .sorted.bam)_assembly_depth.txt && \
    echo "Done. Removing the temporary file" >> ${snakemake_log}
    rm $TMPDIR/$(basename $depth .sorted.bam).txt
done

rm -rf $TMPDIR
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo $runtime