============================
Full installation validation
============================

The fast CI smoke test verifies that ChloroScan installs and that Snakemake
constructs the complete DAG. It does not download the large databases, install
all Conda environments, or execute the complete workflow.

For a complete run, use Linux with sufficient storage (the extracted CAT
database alone requires approximately 85 GB) and memory. Python 3.10 and Conda
are recommended. The CORGI environment defaults to CPU compatibility; set
CHLOROSCAN_INSTALL_CUDA=1 before provisioning rule environments to opt into
CUDA 12.4 packages, only on compatible hosts.

Download and extract the CAT database first. Locate its db and tax directories;
the latter must contain names.dmp. Clone the matching ChloroScan repository,
activate an environment containing ChloroScan and run from the repository root:

.. code-block:: bash

   bash scripts/validate_full_install.sh /data/assembly.fasta \
     /data/contig_depth.tsv sample_01 /data/chloroscan_results \
     /data/CAT_database/db /data/CAT_database/tax 16

The script checks the full workflow DAG, creates all Snakemake Conda
environments using the Conda frontend, executes the full workflow and then
verifies that ChloroScan.done was produced in the output directory.

The Binny reference database is downloaded into the workflow output
resources/A2K_database directory, avoiding writes inside the installed
Python package. Refer to the output logs when a rule fails.

In GitHub Actions, the testing workflow provides a manual provision_envs
option. This tests environment creation only. Run the complete pipeline on a
suitable workstation or HPC node with the full CAT database available.
