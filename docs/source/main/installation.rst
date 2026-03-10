============
Installation
============

We are still upgrading ChloroScan for better plastid MAG recovery, for the best outcomes, make sure to install latest release.

To install the workflow, use pip. Firstly create a virtual environment. Ensure that your machine strictly follows the python version restrictions: > 3.9.

Run this to check if your python version is compatible:

.. code-block:: bash

    python3 --version

If your python version is between 3.9 and 4.0, you can proceed to set up the virtual environment and install the chloroscan.
Run the commands below to set the virtual environment and install the chloroscan.

.. code-block:: bash

    python3 -m venv chloroscan_env
    source chloroscan_env/bin/activate

.. code-block:: bash

    pip3 install chloroscan

Note we didn't upload chloroscan to conda.

Besides installation, ChloroScan requires several databases: one for taxonomy classification by CAT and another for inferring plastid MAG qualities by binning.

CAT database
-------------------

Your choices can be either using the CAT-prepared database which contains all proteins in NCBI's nr database. However, it is ram-taking and time-consuming when running CAT, and requires ~300GB of disk.

To assess plastid contig taxonomy, we recommend downloading our curated protein database consisting of all plastid genome-encoded proteins from NCBI and proteins from UniRef90 (size: ~85GB after unzipping).

command to download: 

.. code-block:: bash

   figshare download -o CAT_db.tar.gz 27990278

Binny database and conda environments
-------------------

You don't need to download binny's database manually, ChloroScan will load it once you installed the conda environment for each job.

Run chloroscan with the flag "``--conda-create-envs-only``" added to setup the conda environments.

.. code-block:: bash

    chloroscan run --conda-create-envs-only

Doing so marks the end of preparation of chloroscan for running. 
