============
Installation
============

We are still upgrading ChloroScan for better plastid MAG recovery, for the best outcomes, make sure to install latest release.

1. pip and conda/mamba
-----------------------

To install the workflow, use pip. Firstly create a virtual environment. We recommend Python 3.10 for the current workflow. The Python package metadata permits Python 3.9 and later, but the complete workflow has not been validated across all those versions. Use Python 3.10 for reproducible setup.

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
    source ~/.bashrc

    conda --version
    mamba --version

    which mamba
    which conda

``NOTE``: Use mamba 1.x (not mamba 2.x), current ChloroScan works well with mamba v1. As a reference, we developed ChloroScan under mamba v1.4.2 and conda v23.3.1.
Now you can proceed to set up the virtual environment and install the chloroscan.

Cargo and FragGeneScanRs should also be installed, for coding sequence extraction, so we prepared a bash script to do this: https://github.com/Andyargueasae/chloroscan/blob/main/install_fraggenescanrs.sh.

.. code-block:: bash

    # Download the script from the repository first, then run it:
    bash install_fraggenescanrs.sh

Simply follow its prompts. Now you are ready to install ChloroScan via pip.

.. code-block:: bash

    python3.10 -m venv chloroscan_env
    source chloroscan_env/bin/activate

.. code-block:: bash

    pip3 install chloroscan

Alternatively, you can create a conda environment with compatible python version, and install chloroscan in it. Make sure that conda/mamba is initiated for bash. 

.. code-block:: bash

    mamba create -n chloroscan_env python=3.10 -y
    conda activate chloroscan_env
    which pip; which python; # check if the pip and python are from the conda environment, 
    #!CAUTION!: default python in your machine shouldn't be shown here, it causes head-scratching issues.
    python -m pip install chloroscan
    # check chloroscan is installed correctly
    chloroscan --help
    chloroscan run --help

Besides installation, ChloroScan requires several databases: one for taxonomy classification by CAT and another for inferring plastid MAG qualities by binning.

2. CAT database
-------------------

Your choices can be either using the CAT-prepared database which contains all proteins in NCBI's nr database. However, it is ram-taking and time-consuming when running CAT, and requires ~300GB of disk.

To assess plastid contig taxonomy, we recommend downloading our curated protein database consisting of all plastid genome-encoded proteins from NCBI and proteins from UniRef90 (size: ~85GB after unzipping). Ensure that you have sufficient disk space for both the archive and extracted files.

command to download: 

.. code-block:: bash

   figshare download -o CAT_db.tar.gz 27990278
   mkdir -p /path/to/your/CAT_database
   tar -xzf CAT_db.tar.gz -C /path/to/your/CAT_database

Identify the extracted CAT database and taxonomy directories (containing files such as ``names.dmp`` in the taxonomy directory). These must be supplied explicitly as ``--cat-database`` and ``--cat-taxonomy``. ChloroScan no longer assumes database paths belonging to the developer.

3. Binny database and conda environments
-----------------------------------------

You don't need to download binny's database manually, ChloroScan will load it once you installed the conda environment for each job.

Run chloroscan with the flag "``--conda-create-envs-only``" added to setup the conda environments. Always add the flag "``--use-conda``" to make sure the workflow uses conda environments. You can specify the prefix of conda environments with "``--conda-prefix``" flag, otherwise the environments will be created in default location. You may always specify the prefix each time you run ChloroScan to avoid rebuilding.

For a typical run, provide ``Inputs-assembly``, ``Inputs-batch-name``, ``outputdir`` and ``Inputs-alignment`` (or a depth file), plus the two required CAT paths ``cat-database`` and ``cat-taxonomy``. (You can use test data to first do the conda environment setup)

.. code-block:: bash

    chloroscan run --use-conda --conda-create-envs-only --conda-prefix="/path/to/your/conda/envs" \
        --Inputs-assembly "/path/to/your/assembly.fasta" --Inputs-batch-name "sample_name" \
        --outputdir "/path/to/your/output_dir" --Inputs-alignment "/path/to/your/bams_directory" \
        --cat-database "/path/to/your/CAT_database/db" \
        --cat-taxonomy "/path/to/your/CAT_database/tax"

The example paths are placeholders: check your extracted database layout before running. This command prepares conda environments only; it does not verify a complete end-to-end run. 
