
===================
Beginner's tutorial
===================

Data Input and configuration
============================

1. example commands.
----------------------------

.. code-block:: bash

   # You have assembly in fasta format and a directory of sorted bam files, you want to run de novo analysis.
   ASSEMBLY="path/to/your/input-contigs.fasta"
   ALIGNMENT="path/to/your/bam_directory"
   BATCH_NAME="your_batch_name"
   OUTPUT="path/to/your/chloroscan_output"
   BINNY_OUTPUTDIR="path/to/your/binny_output_directory"

   chloroscan run --Inputs-assembly=$ASSEMBLY \
      --Inputs-alignment=$ALIGNMENT \
      --Inputs-batch-name=$BATCH_NAME \
      --outputdir=$OUTPUT \
      --use-conda \
      --conda-prefix="conda_env" \
      --conda-frontend="mamba" \
      --cores=12 \
      --verbose \
      --corgi-pthreshold=0.5 \
      --binning-clustering-hdbscan-min-sample-range=1,5,10 \
      --binning-bin-quality-purity=90 \
      --binning-outputdir=$BINNY_OUTPUTDIR \
      --corgi-min-length=1000 \
      --force \
      --cat-database "path/to/your/CAT_db/db" \
      --cat-taxonomy "path/to/your/CAT_db/tax" \
      --binning-bin-quality-min-completeness=50

   # You have assembly in fasta format and a depth profile table, you want to rerun your analysis.
   DEPTH="path/to/your/depth_profile.tsv"
   chloroscan run --Inputs-assembly=$ASSEMBLY \
      --Inputs-depth-profile=$DEPTH \
      --Inputs-batch-name=$BATCH_NAME \
      --outputdir=$OUTPUT \
      --use-conda \
      --conda-prefix="conda_env" \
      --conda-frontend="mamba" \
      --cores=12 \
      --verbose \
      --corgi-pthreshold=0.5 \
      --binning-clustering-hdbscan-min-sample-range=1,5,10
      --binning-bin-quality-purity=90 \
      --binning-outputdir=$BINNY_OUTPUTDIR \
      --corgi-min-length=1000 \
      --force \
      --cat-database "path/to/your/CAT_db/db" \
      --cat-taxonomy "path/to/your/CAT_db/tax" \
      --binning-bin-quality-min-completeness=50


2. Cooperate with Anvi'o to make sense of your database
-------------------------------------------------------

Anvi'o is a powerful tool for visualizing and analyzing metagenomic data. After running ChloroScan, you can use Anvi'o to visualize the results and gain insights into your data.
Potentially it can help you to refine the plastid MAGs. 

The inputs required for Anvi'o are:

- ``reformatted fasta-format assembly``: Filtered contigs from contig classification module, anvi-script-reformat-fasta can take the assembly and output this file.

- ``bam files``: sorted.bam files that contains all contigs in the filtered assembly. Original bam files will fail anvi'o's workflow, so users can use bowtie2/minimap2 to map their raw reads back to the filtered assembly and generate the bam files.

- ``binning collection (tsv)``: The binning mapping file, with the first column being the contig id and the seconda column being the bin id.

After doing these, you can put these data to anvi'o for further analysis. 