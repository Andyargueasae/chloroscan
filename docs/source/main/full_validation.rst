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

Troubleshooting a first run
---------------------------

* **CLI rejects an input flag:** inspect ``chloroscan run --help``; use ``--Inputs-depth-txt`` for an existing depth table, or ``--Inputs-alignment`` for a BAM directory.
* **CAT paths are rejected or rules cannot find taxonomy:** inspect the extracted database hierarchy. ``--cat-database`` must point to the CAT ``db`` directory, and ``--cat-taxonomy`` to a directory containing ``names.dmp``.
* **Conda creation fails:** run with ``--use-conda --conda-frontend conda``; verify channels, disk space, network access and cache permissions. Historic Conda/Mamba downgrades are not required.
* **CUDA/PyTorch installation fails:** leave ``CHLOROSCAN_INSTALL_CUDA`` unset for the default CPU-compatible installation. Only opt into CUDA 12.4 on compatible hardware.
* **A biological workflow rule fails:** check ``OUTPUT_DIR/logging_info`` and Snakemake's error text. Ensure contig identifiers match between FASTA and depth/BAM inputs and that the output directory is writable.

A successful environment-provisioning CI run is **not** evidence of a successful end-to-end biological analysis. Validate at least one complete small dataset on a suitable Linux workstation or HPC system before describing the release as fully tested.
