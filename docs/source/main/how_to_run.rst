
======================
Running ChloroScan
======================

ChloroScan is specialized for plastid genome bin recovery, particularly ``photosynthetic microalgae and macroalgae and complex plastids derived from them, such as diatoms and haptophytes``. The primary inputs it wants are contigs and their depth profiles.
For most cases the default settings will be working well, but you can also configure the parameters for each step. Here we will briefly introduce the inputs and outputs of each step, and the parameters you can configure for each step.

1. Are there any plastids in your data?
======================================

Photosynthetic algae/protists commonly dwell in environments like ocean waters, freshwater, lichen and humid soils. Anthropogenic environments, underground environments and sediments are less likely to have them in high abundances. So if your data come from these environments, ChloroScan may not recover any plastid MAGs.

Meanwhile, if your data contain too fragmented plastid contigs with low coverage (commonly < 5.0 x), ChloroScan may still miss them. Hopefully our newer versions could resolve these issues. 

2. A common minimal command
==========================

You've got your contigs, mapping files in bam and configuration files. You now can run the whole workflow in one command:

.. code-block:: bash

   chloroscan run --Inputs-assembly input_contigs.fasta --Inputs-alignment PATH/to/bams \
      --Inputs-batch-name "my_batch" --outputdir Path/to/output --use-conda --cores=12 \
      --cat-database PATH/to/CAT_db/db --cat-taxonomy path/to/CAT_db/tax

Sometimes you may also wish to get a tabular represented depth profile. It looks like:

.. code-block::

   S0C861	3.52491757
   S0C1664	2.73124830
   S0C2713	12.64139886
   S0C3242	2.51473363
   S0C8106	23.82718202
   S0C8631	2.69335600
   S0C9609	2.49439900
   ...

The first column is the contig id, and the second column is the average depth of contig. We also accept this format with the command changed into:

.. code-block:: bash

   DEPTH_PROFILE=path/to/depth_profile.tsv
   chloroscan run --Inputs-assembly input_contigs.fasta --Inputs-depth-profile $DEPTH_PROFILE \
      --Inputs-batch-name "my_batch" --outputdir Path/to/output --use-conda --cores=12 \
      --cat-database PATH/to/CAT_db/db --cat-taxonomy path/to/CAT_db/tax

For more details about the inputs, please check the :ref:`inputs_and_outputs` section.

If you want to run with our test data, you can use the commands shown in :ref:`README` to download the test data:
  
.. code-block:: bash
   
   figshare download -o simulated_metagenomes.tar.gz 28748540


3. Explanations to arguments in commands
========================================

The whole command space of ChloroScan is shown below:

.. code-block:: 

    Usage: chloroscan run [OPTIONS]                                                            
                                                                                                
    Run the workflow.                                                                          
    All unrecognized arguments are passed onto Snakemake.                                      
                                                                                                
    ╭─ Options ────────────────────────────────────────────────────────────────────────────────╮
    │ --config                   FILE     Path to snakemake config file. Overrides existing    │
    │                                     workflow configuration.                              │
    │                                     [default: None]                                      │
    │ --resource        -r       PATH     Additional resources to copy from workflow directory │
    │                                     at run time.                                         │
    │ --profile         -p       TEXT     Name of profile to use for configuring Snakemake.    │
    │                                     [default: None]                                      │
    │ --dry             -n                Do not execute anything, and display what would be   │
    │                                     done.                                                │
    │ --lock            -l                Lock the working directory.                          │
    │ --dag             -d       PATH     Save directed acyclic graph to file. Must end in     │
    │                                     .pdf, .png or .svg                                   │
    │                                     [default: None]                                      │
    │ --cores           -c       INTEGER  Set the number of cores to use. If None will use all │
    │                                     cores.                                               │
    │                                     [default: None]                                      │
    │ --no-conda                          Do not use conda environments.                       │
    │ --keep-resources                    Keep resources after pipeline completes.             │
    │ --keep-snakemake                    Keep .snakemake folder after pipeline completes.     │
    │ --verbose         -v                Run workflow in verbose mode.                        │
    │ --help-snakemake  -hs               Print the snakemake help and exit.                   │
    │ --help            -h                Show this message and exit.                          │
    ╰──────────────────────────────────────────────────────────────────────────────────────────╯
    ╭─ Workflow Configuration ─────────────────────────────────────────────────────────────────╮
    │ --Inputs-assembly        -a                            PATH     Path to fasta format     │
    │                                                                 assembly of contigs from │
    │                                                                 all sorts of organisms.  │
    │                                                                 [default: None]          │
    │ --Inputs-depth-txt       -d                            PATH     Path to a tab-separated  │
    │                                                                 text storing abundance   │
    │                                                                 of each contig in the    │
    │                                                                 sample.                  │
    │                                                                 [default: None]          │
    │ --Inputs-alignment       -l                            PATH     Path to the folder       │
    │                                                                 containing alignment     │
    │                                                                 files of the contigs.    │
    │                                                                 [default: None]          │
    │ --Inputs-batch-name      -b                            TEXT     Name of the batch.       │
    │                                                                 [default: None]          │
    │ --outputdir              -o                            PATH     Path to the output       │
    │                                                                 directory of the         │
    │                                                                 workflow.                │
    │                                                                 [default: None]          │
    │ --tmpdir                 -t                            PATH     Path to the temporary    │
    │                                                                 directory of the         │
    │                                                                 workflow.                │
    │                                                                 [default: tmp]           │
    │ --binning-universal-le…                                INTEGER  Length cutoff for        │
    │                                                                 universal binning.       │
    │                                                                 [default: 1500]          │
    │ --binning-snakemake-env                                TEXT     Customized snakemake     │
    │                                                                 environment for binny to │
    │                                                                 run.                     │
    │                                                                 [default: None]          │
    │ --binning-mantis-env                                   TEXT     Customized Mantis        │
    │                                                                 virtual environment to   │
    │                                                                 have mantis_pfa          │
    │                                                                 installed, annotating    │
    │                                                                 genes.                   │
    │                                                                 [default: None]          │
    │ --binning-outputdir      -o                            PATH     Path to the output       │
    │                                                                 directory of the         │
    │                                                                 binning.                 │
    │                                                                 [default: binny_output]  │
    │ --binning-clustering-e…                                TEXT     Range of epsilon values  │
    │                                                                 for HDBSCAN clustering.  │
    │                                                                 [default: 0.250,0.000]   │
    │ --binning-clustering-h…                                TEXT     Range of min_samples     │
    │                                                                 values for HDBSCAN       │
    │                                                                 clustering, larger value │
    │                                                                 means larger MAGs.       │
    │                                                                 [default: 1,4,7,10]      │
    │ --binning-bin-quality-…                                FLOAT    Starting completeness    │
    │                                                                 for bin quality.         │
    │                                                                 [default: 92.5]          │
    │ --binning-bin-quality-…                                FLOAT    Minimum completeness for │
    │                                                                 bin quality.             │
    │                                                                 [default: 72.5]          │
    │ --binning-bin-quality-…                                FLOAT    Purity for bin quality.  │
    │                                                                 [default: 95]            │
    │ --corgi-min-length                                     INTEGER  Minimum length of        │
    │                                                                 contigs to be processed  │
    │                                                                 by CORGI.                │
    │                                                                 [default: 500]           │
    │ --corgi-save-filter          --no-corgi-save-filter             Save the filtered        │
    │                                                                 contigs by CORGI (Note:  │
    │                                                                 may take long time).     │
    │                                                                 [default:                │
    │                                                                 no-corgi-save-filter]    │
    │ --corgi-batch-size                                     INTEGER  Batch size for CORGI to  │
    │                                                                 process contigs.         │
    │                                                                 [default: 1]             │
    │ --corgi-pthreshold                                     FLOAT    P-value threshold for    │
    │                                                                 CORGI to determine if    │
    │                                                                 the contigs category is  │
    │                                                                 authentically plastidial │
    │                                                                 or something else.       │
    │                                                                 [default: 0.9]           │
    │ --cat-database           -d                            PATH     Path to the database of  │
    │                                                                 chloroplast genomes.     │
    │                                                                 [default:                │
    │                                                                 /home/yuhtong/scratch/a… │
    │ --cat-taxonomy           -t                            PATH     Path to the taxonomy of  │
    │                                                                 the database.            │
    │                                                                 [default:                │
    │                                                                 /home/yuhtong/scratch/a… │
    │ --krona-env                                            TEXT     Path to the Krona        │
    │                                                                 environment.             │
    │                                                                 [default: kronatools]    │
    ╰──────────────────────────────────────────────────────────────────────────────────────────╯
Below lists those arguments for ChloroScan.
 - ``--Inputs-assembly``: Path to fasta format assembly of contigs from all sorts of organisms.
 - ``--Inputs-depth-txt``: Path to a tab-separated text storing abundance of each contig in the sample. The first column is the contig id, and the second column is the average depth of contig. We also accept this format with the command changed into ``--Inputs-depth-profile``.
 - ``--Inputs-alignment``: Path to the folder containing alignment files of the contigs. The alignment files should be in bam format, and named as "sample_name.bam". The sample name will be extracted from the bam file name by removing the ".bam" suffix. The sample name will be used in the downstream analysis and output files.
 - ``--Inputs-batch-name``: Name of the batch. This will be used in the downstream analysis and output files.
 - ``--outputdir``: Path to the output directory of the workflow. The final results will be stored in this directory. The intermediate results will be stored in a subdirectory called "working" under the output directory. The default value is "output".
 - ``--tmpdir``: Path to the temporary directory of the workflow. The default value is "tmp".
 - ``--binning-universal-length-cutoff``: Length cutoff for universal binning. Contigs shorter than this length will be filtered out before binning. The default value is 1500bp.
 - ``--binning-snakemake-env``: Customized snakemake environment for binny to run. If not specified, the default conda environment will be used.
 - ``--binning-mantis-env``: Customized Mantis virtual environment to have mantis_pfa installed, annotating genes. If not specified, the default conda environment will be used.
 - ``--binning-outputdir``: Path to the output directory of the binning. The default value is "binny_output".
 - ``--binning-clustering-epsilon-range``: Range of epsilon values for HDBSCAN clustering. The default value is "0.250,0.000".
 - ``--binning-clustering-hdbscan-min-sample-range``: Range of min_samples values for HDBSCAN clustering, larger value means larger MAGs. The default value is "1,4,7,10".
 - ``--binning-bin-quality-purity``: Minimum purity for bin quality. The default value is 95.
 - ``--binning-bin-quality-starting-completeness``: Starting completeness for bin quality. Binny uses a sliding completeness to filter bins. The default value is 92.5.
 - ``--binning-bin-quality-min-completeness``: Minimum completeness for bin quality. The default value is 72.5.
 - ``--corgi-min-length``: Minimum length of contigs to be processed by CORGI. The default value is 500bp.
 - ``--corgi-save-filter``: Save the filtered contigs by CORGI (Note: may take long time). The default value is no-corgi-save-filter.
 - ``--corgi-batch-size``: Batch size for CORGI to infer contigs' taxonomic labels. The default value is 1.
 - ``--corgi-pthreshold``: P-value threshold for CORGI to determine if the contigs category is authentically plastidial or something else. The default value is 0.9.
 - ``--cat-database``: Path to the database of CAT saving diamond-processed protein sequences. The default value is "PATH/TO/CAT_db/db".
 - ``--cat-taxonomy``: Path to the taxonomy labels of the CAT database. The default value is "PATH/TO/CAT_db/tax".
 - ``--krona-env``: Path to the Krona environment. The default value is "kronatools". **Now commonly we don't need to set up this**.