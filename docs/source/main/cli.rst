===============================
Command Line Interface
===============================

ChloroScan uses `Snk <https://snk.wytamma.com/snk-cli/>`_ to create a command-line utility. Here are the available commands:

* :ref:`main/cli:run` - Run the ChloroScan pipeline.
* :ref:`main/cli:config` - Show the workflow configuration.
* :ref:`main/cli:env` - Access the workflow conda environments.
* :ref:`main/cli:script` - Access the workflow scripts.
* :ref:`main/cli:info` - Show information about the workflow.
* :ref:`main/cli:profile` - Access the workflow profiles.
* :ref:`main/cli:github` - Launch the ChloroScan GitHub page. 
* :ref:`main/cli:docs` - Launch the ChloroScan documentation.     

To find out more about a command, run:

.. code-block:: bash

    chloroscan --help

run
=========

The main command to run the ChloroScan pipeline is:

.. code-block:: bash

    chloroscan run

The run command will pass all unrecognized arguments to Snakemake. 
That means that if you want to use any of the Snakemake options, you can pass them to the run command e.g. ``chloroscan run --use-singularity``. 
To see all the options available in Snakemake, you can use the ``--help-snakemake`` (``-hs``) flag. 

You can use the ``--config`` flag to specify a configuration file to override the existing workflow configuration. 
This is the same as using the ``--configfile`` flag in Snakemake.

.. code-block:: bash

    chloroscan run --config config.yaml

You can use the ``--profile`` flag to specify a profile to use for configuring Snakemake. You can specify the profile by name. The profile must be located in the profiles directory of the workflow.

.. code-block:: bash

    chloroscan run --profile slurm


You can use the ``--dag`` flag to save the directed acyclic graph to a file. The output file must end in .pdf, .png, or .svg.

.. code-block:: bash

    chloroscan run --dag workflow.pdf

.. note 

    The ``--dag`` flag requires the ``graphviz`` package (``dot``) to be installed.

You can use the ``--dry`` flag to display what would be done without executing anything (this is the same as using the ``--dry-run`` flag in Snakemake).

The ``--lock`` flag is used to lock the working directory (by default, the working directory is not locked e.g. ``--nolock`` is passed to Snakemake).

The ``--cores`` flag is used to set the number of cores to use. If None is specified, all cores will be used by default.

The ``--no-conda`` flag is used to disable the use of conda environments.

The ``--keep-snakemake`` flag is used to keep the .snakemake folder after the pipeline completes. This is useful for debugging purposes. By default, the .snakemake folder is removed after the pipeline completes.

The chloroscan help message is broken into two sections. 
The first section lists the options available for the run command. 
The second section lists the workflow configuration options (generated from the snakemake configfile).

.. code-block:: bash

    chloroscan run --help

That will output something like:

.. code-block:: 

    Usage: chloroscan run [OPTIONS]                                                                                  
                                                                                                                    
    Run the workflow.                                                                                                
    All unrecognized arguments are passed onto Snakemake.                                                            
                                                                                                                    
    ╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────╮
    │ --config                   FILE     Path to snakemake config file. Overrides existing workflow configuration.  │
    │                                     [default: None]                                                            │
    │ --resource        -r       PATH     Additional resources to copy from workflow directory at run time.          │
    │ --profile         -p       TEXT     Name of profile to use for configuring Snakemake. [default: None]          │
    │ --dry             -n                Do not execute anything, and display what would be done.                   │
    │ --lock            -l                Lock the working directory.                                                │
    │ --dag             -d       PATH     Save directed acyclic graph to file. Must end in .pdf, .png or .svg        │
    │                                     [default: None]                                                            │
    │ --cores           -c       INTEGER  Set the number of cores to use. If None will use all cores.                │
    │                                     [default: None]                                                            │
    │ --no-conda                          Do not use conda environments.                                             │
    │ --keep-resources                    Keep resources after pipeline completes.                                   │
    │ --keep-snakemake                    Keep .snakemake folder after pipeline completes.                           │
    │ --verbose         -v                Run workflow in verbose mode.                                              │
    │ --help-snakemake  -hs               Print the snakemake help and exit.                                         │
    │ --help            -h                Show this message and exit.                                                │
    ╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
    ╭─ Workflow Configuration ───────────────────────────────────────────────────────────────────────────────────────╮
    │ --Inputs-assembly                    -a                            PATH     Path to fasta format assembly of   │
    │                                                                             contigs from all sorts of          │
    │                                                                             organisms.                         │
    │                                                                             [default: None]                    │
    │ --Inputs-depth-txt                   -d                            PATH     Path to a tab-separated text       │
    │                                                                             storing abundance of each contig   │
    │                                                                             in the sample.                     │
    │                                                                             [default: None]                    │
    │ --Inputs-alignment                   -l                            PATH     Path to the folder containing      │
    │                                                                             alignment files of the contigs.    │
    │                                                                             [default: None]                    │
    │ --Inputs-batch-name                  -b                            TEXT     Name of the batch. [default: None] │
    │ --outputdir                          -o                            PATH     Path to the output directory of    │
    │                                                                             the workflow.                      │
    │                                                                             [default: None]                    │
    │ --tmpdir                             -t                            PATH     Path to the temporary directory of │
    │                                                                             the workflow.                      │
    │                                                                             [default: tmp]                     │
    │ --binning-universal-length-cutoff                                  INTEGER  Length cutoff for universal        │
    │                                                                             binning.                           │
    │                                                                             [default: 1500]                    │
    │ --binning-snakemake-env                                            TEXT     Customized snakemake environment   │
    │                                                                             for binny to run.                  │
    │                                                                             [default: None]                    │
    │ --binning-mantis-env                                               TEXT     Customized Mantis virtual          │
    │                                                                             environment to have mantis_pfa     │
    │                                                                             installed, annotating genes.       │
    │                                                                             [default: None]                    │
    │ --binning-outputdir                  -o                            PATH     Path to the output directory of    │
    │                                                                             the binning.                       │
    │                                                                             [default: binny_output]            │
    │ --binning-clustering-epsilon-range                                 TEXT     Range of epsilon values for        │
    │                                                                             HDBSCAN clustering.                │
    │                                                                             [default: 0.250,0.000]             │
    │ --binning-clustering-hdbscan-min-s…                                TEXT     Range of min_samples values for    │
    │                                                                             HDBSCAN clustering, larger value   │
    │                                                                             means larger MAGs.                 │
    │                                                                             [default: 1,4,7,10]                │
    │ --binning-bin-quality-starting-com…                                FLOAT    Starting completeness for bin      │
    │                                                                             quality.                           │
    │                                                                             [default: 92.5]                    │
    │ --binning-bin-quality-min-complete…                                FLOAT    Minimum completeness for bin       │
    │                                                                             quality.                           │
    │                                                                             [default: 72.5]                    │
    │ --binning-bin-quality-purity                                       FLOAT    Purity for bin quality.            │
    │                                                                             [default: 95]                      │
    │ --corgi-min-length                                                 INTEGER  Minimum length of contigs to be    │
    │                                                                             processed by CORGI.                │
    │                                                                             [default: 500]                     │
    │ --corgi-save-filter                      --no-corgi-save-filter             Save the filtered contigs by CORGI │
    │                                                                             (Note: may take long time).        │
    │                                                                             [default: no-corgi-save-filter]    │
    │ --corgi-batch-size                                                 INTEGER  Batch size for CORGI to process    │
    │                                                                             contigs.                           │
    │                                                                             [default: 1]                       │
    │ --corgi-pthreshold                                                 FLOAT    P-value threshold for CORGI to     │
    │                                                                             determine if the contigs category  │
    │                                                                             is authentically plastidial or     │
    │                                                                             something else.                    │
    │                                                                             [default: 0.9]                     │
    │ --cat-database                       -d                            PATH     Path to the database of            │
    │                                                                             chloroplast genomes.               │
    │                                                                             [default:                          │
    │                                                                             /home/yuhtong/scratch/andy/uniref… │
    │ --cat-taxonomy                       -t                            PATH     Path to the taxonomy of the        │
    │                                                                             database.                          │
    │                                                                             [default:                          │
    │                                                                             /home/yuhtong/scratch/andy/uniref… │
    │ --krona-env                                                        TEXT     Path to the Krona environment.     │
    │                                                                             [default: kronatools]              │
    ╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

config
=========

The config subcommand will display the workflow configuration file contents. 
You can use the ``--pretty`` (``-p``) flag to display the configuration in a more readable format.

.. code-block:: bash

    chloroscan config

You can redirect the output to a file to save the configuration.

.. code-block:: bash

    chloroscan config > config.yaml

You can then edit the configuration file and use it to run the workflow.

.. code-block:: bash

    chloroscan run --config config.yaml

env
=========

The env subcommand in the workflow tool allows you to access and manage the conda environments used within the workflow. This guide provides an overview of the available options and commands for working with workflow environments.

env list
--------

List the environments in the workflow.

.. code-block:: bash

    chloroscan env list [OPTIONS]

env activate
-------------

This command activates the specified conda environment within the workflow.

.. code-block:: bash

    chloroscan env activate [OPTIONS] ENV_NAME

``ENV_NAME``: Name of the environment to activate.

env create
----------

This command creates all the conda environments specified in the envs dir. Individual conda envs can be create with workflow env create ``ENV_NAME``.

Snakemake workflows that use a lot of conda environments can take a long time to install as each env is created sequentially. Running workflow env create --workers number_of_workers will create all the conda envs in parallel up to the number of workers requested (defaults to 1).

.. code-block:: bash

    chloroscan env create --workers 4  # create up to 4 conda envs at a time

.. warning::

    Some conda envs may not support parallel creation. If you encounter an error, try reducing the number of workers.


env remove
-----------

This command deletes all the conda environments in the workflow. You can also delete individual environments by specifying the environment name. Use the ``--force`` option to skip the confirmation prompt.

.. code-block:: bash

    chloroscan env remove [OPTIONS] [ENV_NAME...]

``ENV_NAME``: Name of the environment to remove.

env run
--------

The env run command in the workflow tool allows you to run a command within one of the workflow environments.

.. code-block:: bash

    chloroscan env run --env ENV_NAME CMD...

``CMD...``: The command and its arguments to execute within the specified environment.

Make sure to replace ``ENV_NAME`` with the actual name of the desired environment, and ``CMD`` with the command you want to run.

For example:

.. code-block:: bash

    chloroscan env run -e my_environment "python script.py --input input_file.txt --output output_file.txt"

This command runs the ``python script.py --input input_file.txt --output output_file.txt`` command within the ``my_environment`` environment in the workflow. 
Adjust the command and environment name according to your specific use case.

env show
--------

This command displays the contents of the environments configuration file used in the workflow.

.. code-block:: bash

    chloroscan env show [OPTIONS]


info
=========

The info subcommand provides JSON formatted information about the workflow.

.. code-block:: bash

    chloroscan info

.. code-block:: json

    {
        "name": "chloroscan",
        "version": "0.1.1",
        "snakefile": "~/chloroscan/workflow/Snakefile",
        "conda_prefix_dir": "~/chloroscan/.conda",
        "singularity_prefix_dir": "~/chloroscan/.singularity",
        "workflow_dir_path": "~/chloroscan"
    }

profile
=========

The profile subcommand provides several commands to manage the workflow profiles. Profiles are used to define different configurations for the workflow e.g. you can configure how the workflow will run on a HPC. You can read more about profiles in the Snakemake documentation.

.. note

    For snk to be able to access the profiles, the profiles must be located in the profiles directory of the workflow.

profile list
-------------

List the profiles in the workflow.

.. code-block:: bash

    chloroscan profile list

profile show
-------------

The show command will display the contents of a profile.

.. code-block:: bash

    chloroscan profile show slurm

You can save the profiles by piping the output to a file.

.. code-block:: bash

    chloroscan profile show slurm > profile/config.yaml

.. note 

    You must save the profile as config.yaml in a directory so that it can be accessed by the workflow.

Load a profile
--------------

You can load a profile by using the ``--profile`` option in the run command.

.. code-block:: bash

    chloroscan run --profile profile/config.yaml

profile edit
--------------

The edit command will open the profile in the default editor. Changes saved will modify the default profile settings for the installation.

.. code-block:: bash

    chloroscan profile edit slurm

script
================

The script commands allow you to interact with the scripts that come with ChloroScan.

script list
------------

The list command will list all scripts in the workflow.

.. code-block:: bash

    chloroscan script list

script show
-------------

The show command will display the contents of a script.

.. code-block:: bash

    chloroscan script show hello.py

script run
----------------

The run command will run a script. Use the ``--env`` option to specify the environment to run the script in.

.. code-block:: bash

    chloroscan script run --env summary merge_contig_depth.py --help


github
==============

.. code-block:: bash

    chloroscan github

Launches the ChloroScan GitHub page at https://github.com/andyargueasae/chloroscan

docs
==============

.. code-block:: bash

    chloroscan docs

Launches the ChloroScan GitHub page at https://andyargueasae.github.io/chloroscan/

