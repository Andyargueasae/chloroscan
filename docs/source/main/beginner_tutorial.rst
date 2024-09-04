
===================
Beginner's tutorial
===================

Data Input and configuration
============================

The pip version of ChloroScan accepts following arguments:

.. image:: docs/source/_static/images/chloroscan_cli.png

Among these arguments, only assembly path, depth text, alignment folder and output folder are required for inputs.  

An example command for running ChloroScan if depth text file is available:

.. code-block:: bash

   chloroscan run --Inputs-assembly-path assembly.fasta --Inputs-depth-txt depth.txt --inputs-batch-name BATCH_NAME --outputdir OUTPUT_DIR

An example command for running ChloroScan if depth text file is not available:

.. code-block:: bash

   chloroscan run --Inputs-assembly-path assembly.fasta --inputs-batch-name BATCH_NAME --outputdir OUTPUT_DIR --alignment-folder ALIGNMENT_FOLDER

Explaining the non-required arguments
=====================================

Some parameters are not required, but changing them can impact the outcomes and running processes, here are some instructions.

1. CORGI settings

CORGI is responsible for contig filtering to extract chloroplast contigs. The default device of CORGI is GPU.

- --corgi-settings-min-length: Minimum length of contigs to be considered for chloroplast contigs. Default is 500. 

- --corgi-settings-pthreshold: P-value threshold for filtering contigs. Default is 0.90. Bigger value means more stringent filtering.

- --corgi-settings-save-filtering: Whether to save the filtered contigs using CORGI. Default is False.

- --corgi-settings-batch-size: number of contigs to be processed in one batch. Default is 1. Smaller batch size can reduce memory usage and elevate speed. 

2. Binning module settings.

- --binning-universal-length-cutoff: The length cutoff for all contigs that enters selection regardless of presenting marker genes. Default is 1500.

- --binning-clustering-epsilon_range: The range of epsilon for HDBSCAN clustering. Default is "0.250,0.000". Larger range takes more time to finish the clustering.

- --binning-clustering-hdbscan-min_samples: The minimum number of samples in a cluster. Default is 1,4,7,10. Picking smaller values can lead to more clusters but cluster contents may be spurious.

- --binning-bin-quality-min_completeness: Minimum completeness for MAGs to enter selection. Default is 72.5.

- --binning-bin-quality-starting-completeness: Starting completeness for MAGs to enter selection. Default is 92.5, note that the range between starting and min completeness affects the quality selections of bins.

- --binning-bin-quality-purity: Purity cutoff for MAGs to enter selection. Default is 95.

3. Homology-based taxonomy prediction by CAT.

- CAT-database: The database used for taxonomy prediction. Default is nrDB, but now there are other databases available such as Uniref90_algaProt.

- CAT-taxonomy: The taxonomy dump file for CAT to use.
