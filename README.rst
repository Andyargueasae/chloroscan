
==================================================================
ChloroScan: A metagenomic workflow to recover chloroplast genomes.
==================================================================


Author: Yuhao Tong (The University of Melbourne)

.. image:: https://github.com/Andyargueasae/chloroscan/badge.svg
   :target: https://github.com/Andyargueasae/chloroscan/badge.svg
   :alt: example workflow

.. image:: docs/source/_static/images/Flow_ChloroScan.drawio.png

A snakemake workflow for MMA metagenomics for recovering chloroplast genomes.
=============================================================================

Project implementation commence date: August 8th, 2023
This workflow will try to cover all essential steps in order to get down the genomics-based metagenomic analysis in recovering algal plastidial genomes.
Before starting the new jobs:

#. Set-up the binny workflow within ChloroScan working directory. 
#. Set-up the CAT-taxonomy identification database, for details please see this: https://tbb.bio.uu.nl/bastiaan/CAT_prepare/
#. Make sure the FragGeneScanRs is added to your path.

All steps above have been covered by running autoInit.sh.

Main dependencies of the workflow:

#. bio-corgi;
#. binny (Customized for ChloroScan);
#. biopython;
#. CAT/BAT (20231120 version requires installation from GitHub, meanwhile the nr database has been updated);
#. diamond;
#. FragGeneScanRs (via cargo and rustc);
#. gffread;
#. numpy;
#. pandas;
#. prodigal

Currently ChloroScan is only available via from-resource installation.

Once the workflow gets finished, the recovered MAGs will be passed to the snakemake workflow Orthoflow (https://rbturnbull.github.io/orthoflow/main/installation.html), to conduct the phylogenetic analysis and see if there any new insights brought by these MAGs.

To install ChloroScan from resource:

.. code-block:: bash
   
   git clone https://github.com/Andyargueasae/chloroscan.git
   cd chloroscan
   poetry install
   poetry shell

Future Update:
1. Comprehensive algae lineage-specific marker gene database.

For more information, please visit the wiki website for details: https://andyargueasae.github.io/chloroscan/. 

