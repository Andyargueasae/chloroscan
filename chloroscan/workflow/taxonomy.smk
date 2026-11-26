rule CAT_taxonomy_identification_plus_annotation:
    input:
        flag=rules.check_binny_output.output.check_flag,
        binny_output=rules.binny_run.output

    params:
        prefix="out.CAT",
        database=CAT_DB_DIR,
        taxonomy=CAT_TAXON_DIR,
        path_to_contigs="intermediary/assembly.formatted.fa"
    output:
        CAT_output=directory(OUTPUT_DIR/"working/CAT")
        # contig_classification=OUTPUT_DIR/"working/CAT/out.CAT.contig2classification.txt"
    log:
        OUTPUT_DIR/"logging_info/CAT.log"
    benchmark:
        OUTPUT_DIR/"logging_info/CAT_benchmark.tsv"
    conda:
        "envs/CAT_env.yml"
    shell:
        """
        set -e
        mkdir -p {output}

        # Check to see if the binny directory is empty or not
        find {OUTPUT_DIR}/working/binny -mindepth 1 -maxdepth 1 >> {log}
        if [ -f "{OUTPUT_DIR}/working/binny/{params.path_to_contigs}" ] ; then
            CAT contigs -c {OUTPUT_DIR}/working/binny/{params.path_to_contigs} -d {params.database:q} -t {params.taxonomy:q} -o {output}/{params.prefix} --force >> {log}
        fi
        """

rule documenting_binny_results:
    input:
        contig_level_annotation = OUTPUT_DIR/"working/CAT"
    params:
        bins = Path(OUTPUT_DIR/"working/binny/bins"),
        assembly_files = OUTPUT_DIR/"working/binny/intermediary/assembly.formatted.fa",
        assembly_depth = OUTPUT_DIR/"working/binny/intermediary/assembly.contig_depth.txt",
        marker_gene_gff = OUTPUT_DIR/"working/binny/intermediary/annotation_CDS_RNA_hmms_checkm.gff",
        names_dump=TAXA,
        MMA_summary=OUTPUT_DIR/"MMA.summary.txt",
        ANNOT_FILE="out.CAT.contig2classification.txt",
        script_dir = SCR_DIR
    output:
        OUTPUT_DIR/"working/summary/cross_ref.tsv",
    benchmark:
        OUTPUT_DIR/"logging_info/summarize_binny.tsv"
    conda:
        "envs/documenting_binny.yml"
    message:
        "Run python script summarize_binny_results.py to record each bin's contig metadata."
    shell:
        """
        $CONDA_PREFIX/bin/python {params.script_dir}/summarize_binny_results.py --assemblyfasta {params.assembly_files} \
            --assemblydepth {params.assembly_depth} --assembly_annot {params.marker_gene_gff} \
            --taxonomy_names {params.names_dump} --contig_level_annotation {input.contig_level_annotation}/{params.ANNOT_FILE} \
            --assembly_bin_dir {params.bins} --MMA_SUMMARY {params.MMA_summary} --output {output}
        """

rule refine_bins:
    input:
        cross_ref = OUTPUT_DIR/"working/summary/cross_ref.tsv",
        original_bins = Path(OUTPUT_DIR/"working/binny"),
        CAT_prediction = Path(OUTPUT_DIR/"working/CAT"),
    params:
        script_dir = SCR_DIR
    output:
        refine_dir = directory(OUTPUT_DIR/"working/refined_bins"),
        summary_info = OUTPUT_DIR/"working/refinement_contig_summary.txt"
    conda:
        "envs/refinement.yml"
    message:
        "Identify potential contaminations in each plastid MAG. Users can decide to use the refined or binny-generated bins."
    log:
        OUTPUT_DIR/"logging_info/refinement.log"
    shell:
        """
        $CONDA_PREFIX/bin/python {params.script_dir}/refine_bins.py --input_cross_ref {input.cross_ref} \
            --input_original_bins {input.original_bins} \
            --input_CAT_prediction {input.CAT_prediction} \
            --output_dir_for_refined_bins {output.refine_dir} \
            --output_summary_info {output.summary_info} \
            --log_file {log}
        """

rule visualize_results:
    input:
        OUTPUT_DIR/"working/summary/cross_ref.tsv",
    params:
        batch_name=BATCH_NAME,
        script_dir = SCR_DIR,
        refine_bins_dir=OUTPUT_DIR/"working/binny",
    output:
        directory(OUTPUT_DIR/"working/visualizations")
    conda:
        "envs/visualization.yml"
    log:
        OUTPUT_DIR/"logging_info/visualizations.log"
    shell:
        """
        $CONDA_PREFIX/bin/python {params.script_dir}/visualization.py --input_cross_ref {input} \
            --batch_name {params.batch_name} \
            --refine_bins_dir {params.refine_bins_dir} \
            --output_dir {output} >> {log} 2>&1
        """

rule krona_taxonomy_plot:
    input:
        CAT_pred_whole_data=rules.CAT_taxonomy_identification_plus_annotation.output.CAT_output,
        # summary_info=rules.documenting_binny_results.output,
    params:
        file_to_remodel="out.CAT.contig2classification.txt",
        file_after_remodel="out.CAT.krona.intake.tsv",
        intake_krona="out.CAT.krona.intake2.tsv",
        script_dir = SCR_DIR,
        batch_name=BATCH_NAME
    conda:
        "envs/krona.yml"
    output:
        OUTPUT_DIR/"Krona.html"
    message:
        "Ensure that in your env there must have pandas, click and of course, Kronatools."
    shell:
        """
        set -e

        # add if statement: if there are empty out.CAT file, we exit with touched files: Krona.html.
        if [ ! -s {input}/{params.file_to_remodel} ]; then
            echo "No contig predictions because previous rules did not find any plastid contigs."
            touch {output}
            exit 0
        fi

        python {params.script_dir}/create_taxonomy_profile.py --cat_prediction_txt {input}/{params.file_to_remodel} --output {input}/{params.file_after_remodel}
        sed 1d {input}/{params.file_after_remodel} > {input}/{params.intake_krona}
        ktImportTaxonomy -m 1 -o {output} {input}/{params.intake_krona}
        """
