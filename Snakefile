'''
Before you run this snakemake workflow, please ensure that:
0. All in all: you must install snakemake at the base environment in order to run this workflow.
1. FragGeneScanRs has been in your path, with intact functionality.
2. Binny snakemake workflow works fine, and the initialization of conda environment and database passed the check pass. 
3. Set up the CAT database.
4. figure out what is gonna be your file's name, put a list into it! [either from snakemake itself or by bash script controlling.]
5. Add CAT's conda environment bin to your path, with intact functionality.
'''

import os
import pandas
import sys
import yaml
import urllib.request
from snakemake.utils import min_version
import logging

min_version("6.0")

shell.executable("bash")

# default configuration file.
configfile: "config/config.default.yaml"
# do we also have to use some string extractions to get the name for batch? 
assembly_path = config['Inputs']['assembly_path']
alignment_files = config['Inputs']['alignment_files']
OUTPUT_DIR = Path(config['outputdir'])
BATCH_NAME = config['Inputs']['batch_name']

BINNY_OUTPUT_DIR = config['binny_settings']['outputdir_binny']    
# For binny: some global names.
GLOBAL_CONFIG = "config.MMA.yaml"
# A functional binny directory is set-up within MMA.
BINNY_DIR = Path("./binny")
SCR_DIR = Path("./script")
# For CAT taxonomy, some global name.
CAT_DB_DIR = Path(config['CAT_database'])

# create the output directory.
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# create temporary directory.
TMP_DIR = config['tmpdir']
if not os.path.isabs(TMP_DIR):
    TMP_DIR = OUTPUT_DIR/TMP_DIR

TMP_DIR.mkdir(parents=True, exist_ok=True) # what if we do want to prepare data and put it here?

# some parameters.

# some functions.

# Now rules are introduced. They are placed into the single snakefile in order to increase efficiency for debugging.
# localrules: ALL 
rule all:
    input:
        # directory(OUTPUT_DIR/"working/binny"),
        # can adjust the input of rule all to choose which rules can be executed. 
        # OUTPUT_DIR/"working/corgi/plastid.fasta"
        # directory(OUTPUT_DIR/"working/binny/cross_ref.xlsx")
        # directory(OUTPUT_DIR/"working/CAT")
        # directory(OUTPUT_DIR/"working/FragGeneScanRs")
        directory(OUTPUT_DIR/"working/gffread")

# The most important two rules are: corgi and binny.
rule corgi_prediction:
    input: 
        seqs = assembly_path 
    params:
        batch_size = config['corgi_settings']['batch_size'],
        min_length = config['corgi_settings']['minlen']
    conda:
        "corgi-env"
    output:
        corgi_out = directory(OUTPUT_DIR/"working/corgi"),
        plastid_contigs = protected(OUTPUT_DIR/"working/corgi/plastid.fasta")
        # here may change the output to a directory for file writing!
    message:
        "Using the machine learning algorithm 'CORGI' to predict contig identity. Don't rise batch size too high."
    resources:
        mem_mb=10000
    threads:
        12
    log:
        OUTPUT_DIR/"logging_info/corgi.log"
    benchmark:
        OUTPUT_DIR/"logging_info/corgi_benchmarking.tsv"
    shell:
        "corgi --file {input.seqs} --output-dir {output.corgi_out} --batch-size {params.batch_size} --min-length {params.min_length} 2> {log}"
        

rule binny_workflow:
    input:
        plastid_contigs = rules.corgi_prediction.output.plastid_contigs,
        # commonly alignment is in the same 
        bamfiles = config['Inputs']['alignment_files']
    params:
        binny_dir = "./binny",
        # outputdir is a must-checked one!
        universal_length_cutoff = config['binny_settings']['universal_length_cutoff'],
        hdbscan_epsilon_range = config['binny_settings']['clustering']['epsilon_range'],
        hdbscan_min_samples_range = config['binny_settings']['clustering']['hdbscan_min_samples_range'],
        quality_min_completeness = config['binny_settings']['bin_quality']['min_completeness'],
        quality_start_completeness = config['binny_settings']['bin_quality']['start_completeness'],
        quality_min_purity = config['binny_settings']['bin_quality']['purity'],
        outputdir_within_binny = config['binny_settings']['outputdir_binny'], # within binny.
         #this place shall be the output, relating to binny!
    output:
        # os.path.join(outputdir, "/working/binny/{params.outputdir_bins_within_binny}/binny.done")
        dir_binny = protected(directory(OUTPUT_DIR/"working/binny")),
        # dir_output = directory(OUTPUT_DIR/"working/binny/bins")
    message:
        "Use the binning algorithm binny to find all high-coverage bins with completeness>=65% and purity>=90%"
    shell:
        """
        bash ./scripts/create_yaml.sh -a {input.plastid_contigs} -l {input.bamfiles} -o {params.outputdir_within_binny}\
           -e {params.hdbscan_epsilon_range} -m {params.hdbscan_min_samples_range} -c {params.quality_min_completeness}\
           -s {params.quality_start_completeness} -p {params.quality_min_purity} -u {params.universal_length_cutoff}
        
        # Now the problem could only be here.
        {params.binny_dir}/binny -l -n "MMA_plastids" -t 16 -r {params.binny_dir}/config/{GLOBAL_CONFIG}

        # must remove the config.yaml created at the end, so we need it to be changed with name.
        # rm {params.binny_dir}/config/{GLOBAL_CONFIG} # finally you have to remove it, you don't want to cause confusion.
        echo 'The binny workflow is done, bins are within the binny folder, that may need to move out.'
        mv {BINNY_OUTPUT_DIR} {output.dir_binny}
        echo 'The binny's output has been moved to the default output directory.'
        """

# rule documenting_binny_results:
#     input: 
#         directory(OUTPUT_DIR/"working/binny/bins")
#     params: 
#         assembly_files = OUTPUT_DIR/"working/binny/intermediary/assembly.formatted.fa",
#         assembly_depth = OUTPUT_DIR/"working/binny/intermediary/assembly.contig_depth.txt",
#         marker_gene_gff = OUTPUT_DIR/"working/binny/intermediary/annotation_CDS_RNA_hmms_checkm.gff"
#     output:
#         OUTPUT_DIR/"working/binny/cross_ref.xlsx"
#     benchmark:
#         OUTPUT_DIR/"logging_info/summarize_binny.tsv"
#     threads:
#         5
#     conda:
#         "base"
#     message:
#         "Run python script summarize_binny_results.py to record each bin's contig information."
#     script:
#         "./scripts/summarize_binny_results.py"


rule CAT_taxonomy_identification_plus_annotation:
    input: 
        bins = directory(OUTPUT_DIR/"working/binny/bins")
    params:
        prefix = config['Inputs']['batch_name']
    output:
        directory(OUTPUT_DIR/"working/CAT"),
    log:
        OUTPUT_DIR/"logging_info/CAT.log"
    benchmark:
        OUTPUT_DIR/"logging_info/CAT_benchmark.tsv"
    threads:
        7
    resources:
        mem_mb=15000
    shell:
        """
        CAT bins --bin_folder {input.bins} -d {CAT_DB_DIR}/2021-01-07_CAT_database -t {CAT_DB_DIR}/2021-01-07_taxonomy -s fasta -o {params.prefix} 
        mkdir -p {output}
        mv ./{params.prefix}.* {output}
        """


rule predict_genes_using_FragGeneScanRs:
    input:
        bins = directory(OUTPUT_DIR/"working/binny/bins")
    params:
        BATCH_NAME
    output:
        # we finally stores all bin's ORF prediction gff files here.
        directory(OUTPUT_DIR/"working/FragGeneScanRs")
    message:
        "Use the fast ORF-predictor FragGeneScan Rust to predict gene contents of the MAGs."
    script:
        "scripts/FragGeneScanRs.sh"


rule gffread_extract_CDS:
    input:
        # we may only want those gff files.
        rules.predict_genes_using_FragGeneScanRs.output
    output:
        directory(OUTPUT_DIR/"working/gffread")
    log:
        OUTPUT_DIR/"logging_info/gffread_CDS.log"
    threads:
        1
    message:
        "Use the gff utility gffread to identify and extract those coding sequences and write the cds file, these will be the input of orthoflow."
    script:
        "scripts/gffread.sh"

# Jobs testing for each of the rule is over, congradulations, now what you should do:
# 1. Set up one complete workflow on a separate virtual machine in order to test its reproducibility.
# 2. test different datasets, especially those larger ones.
# 3. pull those essential files and scripts to github repositories.