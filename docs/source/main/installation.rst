============
Installation
============

Chloroscan can be accessed via github (source code), pypi or conda.

GitHub installation
===================

To install with github, try following (**recommended**):

.. code-block:: bash

    git clone https://github.com/Andyargueasae/chloroscan.git


pip installation
================

Pypi installation of chloroscan is simpler (Developing, do not try now):

.. code-block:: bash

    pip install chloroscan

conda installation
==================

Conda installation (Developing, do not try now):

.. code-block:: bash

    conda install chloroscan

or using mamba with faster speed (Developing, do not try now):

.. code-block:: bash

    mamba install chloroscan

==================================
Key Components' installation
==================================
ChloroScan workflow will be run in snakemake backbone, but before running jobs, make sure to follow the steps below to set it up first. Otherwise workflow lacks its key features to run codes and programs.

I. Individual setups to install and fine-tune all components
=========================================================
a. setting up database for CAT taxonomy assignment. 
CAT is the tool to assign taxon for each contigs inside the dataset you use. Here I recommend using the newest nrDB database.
Alternatively, if you have your preferred versions that are smaller or are at-place, ignore this part and simply specify your database in config file.
Note: Creating a soft link is recommended.

.. code-block:: bash

    # I recommend that you don't install this in ChloroScan. 
    cd PATH/YOU/WANT/TO/SET/YOUR/DATABASE
    wget tbb.bio.uu.nl/tina/CAT_pack_prepare/20231120_CAT_nr.tar.gz
    tar -xvzf 20231120_CAT_nr.tar.gz
    rm 20231120_CAT_nr.tar.gz
    
    # Make a directory inside your ChloroScan working directory. 
    cd PATH/TO/ChloroScan
    mkdir -p databases
    cd databases
    ln -s PATH/OF/YOUR/CAT/DATABASE CAT_db

Voilla. 

b. Set up curl-based packages and add it to $PATH.
In our workflow, the annotator **FragGeneScanRs** could extract cds in fragmented sequences. But it is installable only via cargo. Thus, we need to install it and add it to $PATH, making it available regardless path of your server.
Firstly, install rustc and rustup so that cargo -- the tool to install FragGeneScanRs can be loaded:

.. code-block:: bash
    
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Then, add cargo's path to $PATH.

.. code-block:: bash
    
    . "$HOME/.cargo/env"
    cargo install fraggenescanrs

c. Setting up Krona virtual environment.
Krona is responsible for creating a krona plot that visualizes the abundance of each taxon inside your metagenome dataset. Empirically, it is recommended to create it yourself.

.. code-block:: bash
    
    # Make sure you have conda/mamba installed.
    mamba create -n kronatools -c bioconda -c conda-forge -y krona==2.8.1 pandas numpy

After this, you can specify the conda environment to config file, ChloroScan will incorporate your virtual environment for that job.

d. Set up conda prefix for workflow.
Finally, after sorting out so much things, we can now ask snakemake to set up the rest of environments for us. Simply run:

.. code-block:: bash

    snakemake -c N_CORES --use-conda --conda-prefix=./conda --snakefile ./chloroscan/workflow/Snakefile --configfile config/ChloroScan.init.yaml --conda-create-envs-only

Hurrah! Now you are good to go. 


II. Run it.
===========
Snakemake command line can execute your jobs, just prepare your **contig assembly in fasta format**, and **contig abundance metrics in bam format**, fill in the config yaml file, and run:

.. code-block:: bash

    snakemake -c N_CORES --configfile PATH/TO/CONFIG.yaml --snakefile chloroscan/workflow/Snakefile --use-conda --conda-prefix YOUR/conda 


ChloroScan is also able to run via snk-cli, but firstly it requires poetry to be present. Just type the code below:

.. code-block:: bash

    cd PATH/TO/chloroscan
    poetry install
    # So that all virtual env packages are ready.

    poetry shell # activate virtualenv.

    chloroscan -h # Test using this command.
    chloroscan run --config=PATH/TO/CONFIGFILE --use-conda --conda-prefix="PATH/TO/CONDA" --cores=N_CORES # Run workflow in this way.

For more advanced information, visit the "beginner's tutorial". 