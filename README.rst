
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

To install the workflow, use pip3:

.. code-block:: bash

    pip3 install chloroscan==0.1.5

Detailed workflow instructions can be found at: https://andyargueasae.github.io/chloroscan/index.html

Machine/OS Requirements
=======================
ChloroScan is only tested on Linux (x86_64), running on IOS system is not recommended.


Configuration databases
=======================
Before running ChloroScan, some packages and datasets need to be installed to run CAT taxonomy prediction properly.

ChloroScan incorporates a marker gene database while running binning, you don't need to do anything, it will be loaded since you build conda environments.

To download our curated Uniref90-algae plastid protein database, use the link: https://doi.org/10.26188/27990278. 

To avoid authentication issues, try the command below to download it:

.. code block:: bash

    wget --referer=https://figshare.unimelb.edu.au --user-agent="Mozilla/5.0" -O "CAT_db.tar.gz" https://figshare.unimelb.edu.au/ndownloader/files/51053993

 - Note: The tar.gz format of CAT database's size is 47GB, and nearly 85GB after unzipped, please ensure you have enough disk storage. 

Sample data to try
==================
To try ChloroScan, I recommend downloading our synthetic metagenome data via the command: 

.. code block:: bash

    wget --referer=https://figshare.unimelb.edu.au --user-agent="Mozilla/5.0" -O "synthetic_sample_II.tar.gz" https://figshare.unimelb.edu.au/ndownloader/files/53499140

Credit
============

ChloroScan is developed by:

.. start-credits

- Yuhao Tong (University of Melbourne)
- `Dr Robert Turnbull <https://findanexpert.unimelb.edu.au/profile/877006-robert-turnbull>`_ 
- `Dr Vanessa Rossetto Marcelino <https://findanexpert.unimelb.edu.au/profile/532755-vanessa-rossetto-marcelino>`_ 
- `A/Prof Heroen Verbruggen <https://hverbruggen.github.io/>`_

.. end-credits

