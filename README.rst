
==================================================================
ChloroScan: A metagenomic workflow to recover chloroplast genomes
==================================================================

.. start-badges

|testing badge| |docs badge|

.. |testing badge| image:: https://github.com/Andyargueasae/chloroscan/actions/workflows/testing.yml/badge.svg
    :target: https://github.com/Andyargueasae/chloroscan/actions

.. |docs badge| image:: https://github.com/Andyargueasae/chloroscan/actions/workflows/docs.yml/badge.svg
    :target: https://Andyargueasae.github.io/chloroscan
    
.. end-badges


.. image:: docs/source/_static/images/new_ChloroScan_workflow.drawio.png

This workflow is designed to recover chloroplast genomes from metagenomic datasets.

Installation
============

To install the workflow, use pip3. The background environment will require Python <4.0, >=3.9 to set up the virtual environment.

.. code-block:: bash

    pip3 install chloroscan==0.1.5

Detailed workflow instructions can be found at: https://andyargueasae.github.io/chloroscan/index.html

Machine/OS Requirements
=======================
ChloroScan is only tested on Linux (x86_64), running on IOS system is not recommended. 
ChloroScan can be installed on servers with hpc clusters and it is recommended to use a GPU to accelerate its running.

``Note``: Through testing, current version of chloroscan cannot support NVIDIA H-100 GPU, due to cuda version incompatibilities. We will work on updating it to allow better performances. 


Configuration databases
=======================
Before running ChloroScan, some packages and datasets need to be installed to run CAT taxonomy prediction properly.

ChloroScan incorporates a marker gene database while running binning, you don't need to do anything, it will be loaded since you build conda environments.

To download our curated Uniref90-algae plastid protein database, use the link: https://doi.org/10.26188/27990278. 

To avoid authentication issues, we recommend using the pyfigshare command-line tool to download. The information of this tool can be found at: https://pypi.org/project/pyfigshare/. Python > 3.0 is required to download it.

Before downloading the files, set up your own figshare account and add an api token to the file ~/.figshare/token.
Then run:

.. code-block:: bash

    figshare download -o CAT_db.tar.gz 27990278


``Note``: The tar.gz format of CAT database's size is 47GB, and nearly 85GB after unzipped, please ensure you have enough disk storage. Meanwhile, the space to setup the conda environment also requires 15 GB of disk.  

Sample data to try
==================
To try ChloroScan, I recommend downloading our synthetic metagenome data via the command: 

.. code-block:: bash

    figshare download -o simulated_metagenomes.tar.gz 28748540

There are also some real metagenome datasets (modified to keep them lightweight) available at: https://figshare.unimelb.edu.au/articles/dataset/ChloroScan_test_data/30218614.

To download:

.. code-block:: bash

    figshare download -o real_test_samples.tar.gz 30218614

Credit
============

ChloroScan is developed by:

.. start-credits

- Yuhao Tong (University of Melbourne)
- `Dr Robert Turnbull <https://findanexpert.unimelb.edu.au/profile/877006-robert-turnbull>`_ 
- `Dr Vanessa Rossetto Marcelino <https://findanexpert.unimelb.edu.au/profile/532755-vanessa-rossetto-marcelino>`_ 
- `A/Prof Heroen Verbruggen <https://hverbruggen.github.io/>`_

.. end-credits

