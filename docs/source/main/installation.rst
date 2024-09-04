
==================================================================
ChloroScan: A metagenomic workflow to recover chloroplast genomes.
==================================================================


Author: Yuhao Tong (The University of Melbourne)

.. image:: https://github.com/Andyargueasae/chloroscan/badge.svg
   :target: https://github.com/Andyargueasae/chloroscan/badge.svg
   :alt: example workflow

.. image:: docs/source/_static/images/Flow_ChloroScan.drawio.png

ChloroScan will try to cover all essential steps in order to get down the genomics-based metagenomic analysis in recovering algal plastidial genomes.

Once the workflow gets finished, the recovered MAGs can be passed to the snakemake workflow Orthoflow (https://rbturnbull.github.io/orthoflow/main/installation.html), to perform phylogenetic analysis and gene annotation.


installation
============

To install ChloroScan from resource, first prepare a virtual environment with snakemake installed:

.. code-block:: bash

   mamba create -n chloroscan-env -y snakemake 

Then you can clone it from resources and create conda environments required.

.. code-block:: bash

   git clone https://github.com/Andyargueasae/chloroscan.git
   cd chloroscan
   
   snakemake -c 20 --configfile config/ChloroScan.init.yaml --snakefile chloroscan/workflow/Snakefile --use-conda --conda-create-envs-only --conda-prefix conda

Alternatively, you can install ChloroScan from PyPI:

.. code-block:: bash

   pip install chloroscan



