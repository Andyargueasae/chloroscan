============
Installation
============

Chloroscan can be accessed via github (source code), pypi or conda.

GitHub installation
===================

To install with github, try following:

.. code-block:: bash

    git clone https://github.com/Andyargueasae/chloroscan.git


pip installation
================

Pypi installation of chloroscan is simpler:

.. code-block:: bash

    pip install chloroscan

conda installation
==================

Conda installation:

.. code-block:: bash

    conda install chloroscan

or using mamba with faster speed:

.. code-block:: bash

    mamba install chloroscan

Other Key Components' installation
==================================
Due to the fact that ChloroScan workflow is composed of different packages and some are also snakemake workflows, this may take a while to collect them and install.
Luckily, the bash script autoInit.sh could do this in one piece.

But before running it, downloading the non-redundant protein database in prior to your data storing directory is recommended: tbb.bio.uu.nl/tina/CAT_pack_prepare/20231120_CAT_nr.tar.gz
If not, this script can handle it as well.

.. code-block:: bash

    bash autoInit.sh -a "ABSOLUTE/PATH/TO/YOUR/CHLOROSCAN/DIRECTORY" -b "PATH/TO/YOUR/HMM/PF/DB/FOR/BINNY" -c "PATH/TO/YOUR/TAXON/MARKER/SET/TSV/FOR/BINNY" -d "PATH/TO/YOUR/DIRECTORY/WHERE/DATABASES/ARE"

Note: this process can be run with obstacles. If so, here are individual commands that could finish setting up the workflow, but with more steps to follow.


Individual setups to install and fine-tune all components
=========================================================
a. setting up binny
Binny is a snakemake-workflow that clusters contigs into bins, which works as the core component here. Any problems in setting up this tool could lead to more troubles in downstream. 
We firstly install it:

.. code-block:: bash

    git clone https://github.com/Andyargueasae/binny_Chloroscan.git
    cd binny_ChloroScan

Because originally binny was only used for bacterial sequences, once initializing binny:

.. code-block:: bash

    ./binny -i -r config/config.init.yaml

We have to change the database it uses and the scoring of the different databases it uses. First, you can see there is a folder called **binny_algal_database**, it stores one hmm profile for 56 algal marker genes as well as a tsv file showing two lines of algal marker genes, working to determine the minimal algal plastid genomes.
We need to add these two files to the **database** folder inside binny_ChloroScan:

.. code-block:: bash
    cd binny_ChloroScan/database/hmms/checkm_pf && rm checkm_pf_filtered.hmm
    cp ../../../../binny_algal_database/checkm_pf_filtered.hmm . && hmmpress -f checkm_pf_filtered.hmm
    # Make sure you have hmmer3 installed to your virtual environment.

Then, we shall also change the chunks, and change the taxon_marker_sets_lineage_sorted.tsv file. 
.. code-block:: bash
    # Clear chunks first.
    rm chunks/* 
    # Copy algal hmm profile to chunks dir, working as the chunk_0.
    cp checkm_pf_filtered.hmm chunks/checkm_pf_filtered_chunk_0.hmm && hmmpress -f chunks/checkm_pf_filtered_chunk_0.hmm
    cd ../../ 
    cp ../binny_algal_database/taxon_marker_sets_lineage_sorted.tsv .

Great. Now the next step is to change the config file: binny_mantis.cfg for annotation tool "mantis" to use these databases. You can simply delete the line: **custom_ref="PATH/TO/binny_ChloroScan/database/hmms/checkm_tf/checkm_tf_filtered.hmm"** and its **weight** used by mantis, as this hmm file still stores bacteria markers that binny orignially retrieved. 
Once you delete these two lines, binny will no longer use checkm_tf to annotate genes, it only annotates algal plastid genes.

Great. Moving to the next step.

b. setting up database for CAT taxonomy assignment. 
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

c. Set up curl-based packages and add it to $PATH.
In our workflow, the annotator **FragGeneScanRs** could extract cds in fragmented sequences. But it is installable only via cargo. Thus, we need to install it and add it to $PATH, making it available regardless path of your server.
Firstly, install rustc and rustup so that cargo -- the tool to install FragGeneScanRs can be loaded:

.. code-block:: bash
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Then, add cargo's path to $PATH.
.. code-block:: bash
    . "$HOME/.cargo/env"
    cargo install fraggenescanrs

d. Setting up Krona virtual environment.
Krona is responsible for creating a krona plot that visualizes the abundance of each taxon inside your metagenome dataset. Empirically, it is recommended to create it yourself.
.. code-block:: bash
    # Make sure you have conda/mamba installed.
    mamba create -n kronatools -c bioconda -c conda-forge -y krona==2.8.1 pandas numpy

After this, you can specify the conda environment to config file, ChloroScan will incorporate your virtual environment for that job.

e. Set up conda prefix for workflow.
Finally, after sorting out so much things, we can now ask snakemake to set up the rest of environments for us. Simply run:
.. code-block:: bash
    snakemake -c N_CORES --use-conda --conda-prefix=./conda --snakefile ./chloroscan/workflow/Snakefile --configfile config/ChloroScan.init.yaml --conda-create-envs-only

Hurrah! Now you are good to go. 

ChloroScan is able to be run via snk-cli, but firstly it requires poetry to be present. Just type the code below:

.. code-block:: bash

    cd PATH/TO/chloroscan
    poetry install
    # So that all virtual env packages are ready.

    poetry shell # activate virtualenv.

    chloroscan -h # Test using this command.
    chloroscan run --config=PATH/TO/CONFIGFILE --use-conda --conda-prefix="PATH/TO/CONDA" --cores=N_CORES # Run workflow in this way.
