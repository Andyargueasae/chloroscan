============
Installation
============

Chloroscan can be accessed via github (source code), pypi or conda.

GitHub installation
===================

To install with github, try following:
.. code_block:: bash
    git clone https://github.com/Andyargueasae/chloroscan.git

pip installation
================

Pypi installation of chloroscan is simpler:
.. code_block:: bash
    pip install chloroscan

conda installation
==================

Conda installation:
.. code_block:: bash
    conda install chloroscan

or using mamba with faster speed:
.. code_block:: bash
    mamba install chloroscan

Other Key Components' installation
==================================
Due to the fact that ChloroScan workflow is composed of different packages and some are also snakemake workflows, this may take a while to collect them and install.
Luckily, the bash script autoInit.sh could do this in one piece.

But before running it, downloading the non-redundant protein database in prior to your data storing directory is recommended: tbb.bio.uu.nl/tina/CAT_pack_prepare/20231120_CAT_nr.tar.gz
If not, this script can handle it as well.
(shall change this to argument-assigning.)

.. code_block:: bash
    bash autoInit.sh -a "ABSOLUTE/PATH/TO/YOUR/CHLOROSCAN/DIRECTORY" -b "PATH/TO/YOUR/HMM/PF/DB/FOR/BINNY" -c "PATH/TO/YOUR/TAXON/MARKER/SET/TSV/FOR/BINNY" -d "PATH/TO/YOUR/DIRECTORY/WHERE/DATABASES/ARE"

Now you are good to go.