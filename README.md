===================
# ChloroScan: A metagenomic workflow to recover chloroplast genomes.
===================
## Author: Yuhao Tong (The University of Melbourne)

![example workflow](https://github.com/Andyargueasae/chloroscan/badge.svg)

The collection of snakemake workflow for MMA metagenomics for recovering chloroplast genomes.

Project implementation commence date: August 8th, 2023

This workflow will try to cover all essential steps in order to get down the genomics-based metagenomic analysis in recovering algal plastidial genomes.

Before starting the new jobs:
1. Set-up the binny workflow, with running ./binny -i config/config.init.yaml
2. Set-up the CAT-taxonomy identification database, for details please see this: https://tbb.bio.uu.nl/bastiaan/CAT_prepare/
3. Make sure the FragGeneScanRs is added to your path.

Main dependencies of the workflow:
1. bio-corgi;
2. binny;
3. biopython;
4. CAT/BAT;
5. diamond;
6. FragGeneScanRs (via cargo and rustc);
7. gffread;
8. numpy;
9. pandas;
10. prodigal

Once the workflow gets finished, the recovered MAGs will be passed to the snakemake workflow Orthoflow (https://rbturnbull.github.io/orthoflow/main/installation.html), to conduct the phylogenetic analysis and see if there any new insights brought by these MAGs.

# The results of Orthoflow: stored in this repository as supermatrix.cds.fa.treefile, with bootstrap support.
