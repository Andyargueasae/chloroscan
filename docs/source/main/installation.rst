============
Installation
============

We are still upgrading ChloroScan for better plastid MAG recovery, for the best outcomes, make sure to install latest release.

1. pip and conda/mamba
-----------------------

To install the workflow, use pip. Firstly create a virtual environment. Ensure that your machine strictly follows the python version restrictions: > 3.9 and < 3.12. As a reference, we developed ChloroScan under Python 3.10.4.

Run this to check if your python version is compatible:

.. code-block:: bash

    python3 --version

Another key component of using ChloroScan is an appropriate conda distribution, so you need to install conda/mamba via downloading miniforge from the link: https://github.com/conda-forge/miniforge/releases/. This will give you minimal installers of conda and mamba.

Once you download the installer, run the command below to install mamba (replace the {OPERATING_SYSTEM} and {ARCHITECTURE} with your own system information, e.g., Miniforge3-Linux-x86_64.sh).:

.. code-block:: bash

    bash Miniforge3-{OPERATING_SYSTEM}-{ARCHITECTURE}.sh

Check whether conda and mamba are intactly installed by running the commands below:

.. code-block:: bash

    # initialize conda and mamba for bash.
    conda init bash
    mamba init bash
    source ~/.bashrc

    conda --version
    mamba --version

    which mamba
    which conda

``NOTE``: Make sure mamba's version isn't > 2.0, current ChloroScan works well with mamba v1. As a reference, we developed ChloroScan under mamba v1.4.2 and conda v23.3.1.
Now you can proceed to set up the virtual environment and install the chloroscan.

Cargo and FragGeneScanRs should also be installed, for coding sequence extraction, so we prepared a bash script to do this: https://github.com/Andyargueasae/chloroscan/blob/release_v0.1.7/install_fraggenescanrs.sh.

.. code-block:: bash

    bash install_fraggenescanrs.sh

Simply follow its prompts. Now you are ready to install ChloroScan via pip.

.. code-block:: bash

    python3 -m venv chloroscan_env
    source chloroscan_env/bin/activate

.. code-block:: bash

    pip3 install chloroscan

Alternatively, you can create a conda environment with compatible python version, and install chloroscan in it. Make sure that conda/mamba is initiated for bash. 

.. code-block:: bash

    mamba env create -n chloroscan_env python=3.10
    conda activate chloroscan_env
    which pip; which python; # check if the pip and python are from the conda environment, 
    #!CAUTION!: default python in your machine shouldn't be shown here, it causes head-scratching issues.
    pip3 install chloroscan
    # check chloroscan is installed correctly
    chloroscan --help
    chloroscan run --help

Besides installation, ChloroScan requires several databases: one for taxonomy classification by CAT and another for inferring plastid MAG qualities by binning.

2. CAT database
-------------------

Your choices can be either using the CAT-prepared database which contains all proteins in NCBI's nr database. However, it is ram-taking and time-consuming when running CAT, and requires ~300GB of disk.

To assess plastid contig taxonomy, we recommend downloading our curated protein database consisting of all plastid genome-encoded proteins from NCBI and proteins from UniRef90 (size: ~85GB after unzipping).

command to download: 

.. code-block:: bash

   figshare download -o CAT_db.tar.gz 27990278

3. Binny database and conda environments
-----------------------------------------

You don't need to download binny's database manually, ChloroScan will load it once you installed the conda environment for each job.

Run chloroscan with the flag "``--conda-create-envs-only``" added to setup the conda environments. Always add the flag "``--use-conda``" to make sure the workflow uses conda environments. You can specify the prefix of conda environments with "``--conda-prefix``" flag, otherwise the environments will be created in default location. You may always specify the prefix each time you run ChloroScan to avoid rebuilding.

Current version of ChloroScan requires four mandatory parameters: ``Inputs-assembly``, ``Inputs-batch-name``, ``outputdir`` and ``Inputs-alignment``. (You can use test data to first do the conda environment setup)

.. code-block:: bash

    chloroscan run --use-conda --conda-create-envs-only --conda-prefix="/path/to/your/conda/envs" \
        --Inputs-assembly "/path/to/your/assembly.fasta" --Inputs-batch-name "sample_name" \
        --outputdir "/path/to/your/output_dir" --Inputs-alignment "/path/to/your/bams_directory"

Doing so marks the end of preparation of chloroscan for running. 
