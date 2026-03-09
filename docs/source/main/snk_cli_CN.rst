=========================================
Snakemake-相关的指令
=========================================

除了用于运行宏基因组分析的命令外，ChloroScan 还使用 `Snk <https://snk.wytamma.com/snk-cli/>`_ 创建命令行工具。以下是可用的命令：

* :ref:`main/cli:run` - 运行 ChloroScan 流程。
* :ref:`main/cli:config` - 显示工作流程配置。
* :ref:`main/cli:env` - 访问工作流程的 conda 环境。
* :ref:`main/cli:script` - 访问工作流程脚本。
* :ref:`main/cli:info` - 显示工作流程信息。
* :ref:`main/cli:profile` - 访问工作流程配置文件。
* :ref:`main/cli:github` - 打开 ChloroScan 的 GitHub 页面。
* :ref:`main/cli:docs` - 打开 ChloroScan 文档。

要了解有关命令的更多信息，请运行：

.. code-block:: bash

    chloroscan --help

=========

运行 ChloroScan 流程的主要命令是：

.. code-block:: bash

    chloroscan run

跑这个 run 命令会将所有未识别的参数传递给 Snakemake。
这意味着如果您想使用任何 Snakemake 选项，可以将它们传递给 run 命令，例如 ``chloroscan run --use-singularity``。
要查看 Snakemake 中可用的所有选项，可以使用 ``--help-snakemake``（``-hs``）标志。

您可以使用 ``--config`` 标志指定配置文件以覆盖现有的工作流程配置。
这与在 Snakemake 中使用 ``--configfile`` 标志相同。

.. code-block:: bash

    chloroscan run --config config.yaml

您可以使用 ``--profile`` 标志指定用于配置 Snakemake 的配置文件。您可以按名称指定配置文件。配置文件必须位于工作流程的 profiles 目录中。

.. code-block:: bash

    chloroscan run --profile slurm


您可以使用 ``--dag`` 标志将有向无环图保存到文件中。输出文件必须以 .pdf、.png 或 .svg 结尾。

.. code-block:: bash

    chloroscan run --dag workflow.pdf

.. note 

    ``--dag`` 标志需要安装 ``graphviz`` 包（``dot``）。

您可以使用 ``--dry`` 标志显示将要执行的操作，而不实际执行任何操作（这与在 Snakemake 中使用 ``--dry-run`` 标志相同）。

``--lock`` 标志用于锁定工作目录（默认情况下，工作目录未锁定，例如 ``--nolock`` 被传递给 Snakemake）。

``--cores`` 标志用于设置要使用的核心数。如果指定为 None，默认情况下将使用所有核心。

``--no-conda`` 标志用于禁用 conda 环境的使用。

``--keep-snakemake`` 标志用于在管道完成后保留 .snakemake 文件夹。这对于调试非常有用。默认情况下，管道完成后会删除 .snakemake 文件夹。
ChloroScan 的帮助信息分为两个部分。
第一部分列出了 run 命令可用的选项。
第二部分列出了工作流程配置选项（由 snakemake 配置文件生成）。

.. code-block:: bash

    chloroscan run --help

这会产生以下输出：

.. code-block:: 

    Usage: chloroscan run [OPTIONS]                                                            
                                                                                                
    Run the workflow.                                                                          
    All unrecognized arguments are passed onto Snakemake.                                      
                                                                                                
    ╭─ Options ────────────────────────────────────────────────────────────────────────────────╮
    │ --config                   FILE     Path to snakemake config file. Overrides existing    │
    │                                     workflow configuration.                              │
    │                                     [default: None]                                      │
    │ --resource        -r       PATH     Additional resources to copy from workflow directory │
    │                                     at run time.                                         │
    │ --profile         -p       TEXT     Name of profile to use for configuring Snakemake.    │
    │                                     [default: None]                                      │
    │ --dry             -n                Do not execute anything, and display what would be   │
    │                                     done.                                                │
    │ --lock            -l                Lock the working directory.                          │
    │ --dag             -d       PATH     Save directed acyclic graph to file. Must end in     │
    │                                     .pdf, .png or .svg                                   │
    │                                     [default: None]                                      │
    │ --cores           -c       INTEGER  Set the number of cores to use. If None will use all │
    │                                     cores.                                               │
    │                                     [default: None]                                      │
    │ --no-conda                          Do not use conda environments.                       │
    │ --keep-resources                    Keep resources after pipeline completes.             │
    │ --keep-snakemake                    Keep .snakemake folder after pipeline completes.     │
    │ --verbose         -v                Run workflow in verbose mode.                        │
    │ --help-snakemake  -hs               Print the snakemake help and exit.                   │
    │ --help            -h                Show this message and exit.                          │
    ╰──────────────────────────────────────────────────────────────────────────────────────────╯
    ╭─ Workflow Configuration ─────────────────────────────────────────────────────────────────╮
    │ --Inputs-assembly        -a                            PATH     Path to fasta format     │
    │                                                                 assembly of contigs from │
    │                                                                 all sorts of organisms.  │
    │                                                                 [default: None]          │
    │ --Inputs-depth-txt       -d                            PATH     Path to a tab-separated  │
    │                                                                 text storing abundance   │
    │                                                                 of each contig in the    │
    │                                                                 sample.                  │
    │                                                                 [default: None]          │
    │ --Inputs-alignment       -l                            PATH     Path to the folder       │
    │                                                                 containing alignment     │
    │                                                                 files of the contigs.    │
    │                                                                 [default: None]          │
    │ --Inputs-batch-name      -b                            TEXT     Name of the batch.       │
    │                                                                 [default: None]          │
    │ --outputdir              -o                            PATH     Path to the output       │
    │                                                                 directory of the         │
    │                                                                 workflow.                │
    │                                                                 [default: None]          │
    │ --tmpdir                 -t                            PATH     Path to the temporary    │
    │                                                                 directory of the         │
    │                                                                 workflow.                │
    │                                                                 [default: tmp]           │
    │ --binning-universal-le…                                INTEGER  Length cutoff for        │
    │                                                                 universal binning.       │
    │                                                                 [default: 1500]          │
    │ --binning-snakemake-env                                TEXT     Customized snakemake     │
    │                                                                 environment for binny to │
    │                                                                 run.                     │
    │                                                                 [default: None]          │
    │ --binning-mantis-env                                   TEXT     Customized Mantis        │
    │                                                                 virtual environment to   │
    │                                                                 have mantis_pfa          │
    │                                                                 installed, annotating    │
    │                                                                 genes.                   │
    │                                                                 [default: None]          │
    │ --binning-outputdir      -o                            PATH     Path to the output       │
    │                                                                 directory of the         │
    │                                                                 binning.                 │
    │                                                                 [default: binny_output]  │
    │ --binning-clustering-e…                                TEXT     Range of epsilon values  │
    │                                                                 for HDBSCAN clustering.  │
    │                                                                 [default: 0.250,0.000]   │
    │ --binning-clustering-h…                                TEXT     Range of min_samples     │
    │                                                                 values for HDBSCAN       │
    │                                                                 clustering, larger value │
    │                                                                 means larger MAGs.       │
    │                                                                 [default: 1,4,7,10]      │
    │ --binning-bin-quality-…                                FLOAT    Starting completeness    │
    │                                                                 for bin quality.         │
    │                                                                 [default: 92.5]          │
    │ --binning-bin-quality-…                                FLOAT    Minimum completeness for │
    │                                                                 bin quality.             │
    │                                                                 [default: 72.5]          │
    │ --binning-bin-quality-…                                FLOAT    Purity for bin quality.  │
    │                                                                 [default: 95]            │
    │ --corgi-min-length                                     INTEGER  Minimum length of        │
    │                                                                 contigs to be processed  │
    │                                                                 by CORGI.                │
    │                                                                 [default: 500]           │
    │ --corgi-save-filter          --no-corgi-save-filter             Save the filtered        │
    │                                                                 contigs by CORGI (Note:  │
    │                                                                 may take long time).     │
    │                                                                 [default:                │
    │                                                                 no-corgi-save-filter]    │
    │ --corgi-batch-size                                     INTEGER  Batch size for CORGI to  │
    │                                                                 process contigs.         │
    │                                                                 [default: 1]             │
    │ --corgi-pthreshold                                     FLOAT    P-value threshold for    │
    │                                                                 CORGI to determine if    │
    │                                                                 the contigs category is  │
    │                                                                 authentically plastidial │
    │                                                                 or something else.       │
    │                                                                 [default: 0.9]           │
    │ --cat-database           -d                            PATH     Path to the database of  │
    │                                                                 chloroplast genomes.     │
    │                                                                 [default:                │
    │                                                                 /home/yuhtong/scratch/a… │
    │ --cat-taxonomy           -t                            PATH     Path to the taxonomy of  │
    │                                                                 the database.            │
    │                                                                 [default:                │
    │                                                                 /home/yuhtong/scratch/a… │
    │ --krona-env                                            TEXT     Path to the Krona        │
    │                                                                 environment.             │
    │                                                                 [default: kronatools]    │
    ╰──────────────────────────────────────────────────────────────────────────────────────────╯

配置/config
=========

config指令可以展示工作流程的配置文件内容，您可以使用 ``--pretty``（``-p``）标志以更易读的格式显示配置。

.. code-block:: bash

    chloroscan config

您可以将输出重定向到文件以保存配置。

.. code-block:: bash

    chloroscan config > config.yaml

然后，您可以编辑配置文件并使用它来运行工作流程。

.. code-block:: bash

    chloroscan run --config config.yaml

env
=========

env 子命令允许您访问和管理工作流程中使用的 conda 环境。本指南提供了有关使用工作流程环境的可用选项和命令的概述。

env list
--------

列出工作流程中的环境。

.. code-block:: bash

    chloroscan env list [OPTIONS]

env activate
-------------

此命令在工作流程中激活指定的 conda 环境。
.. code-block:: bash

    chloroscan env activate [OPTIONS] ENV_NAME

``ENV_NAME``: 要激活的环境的名称。

env create
----------

此命令会创建 envs 目录中指定的所有 conda 环境。可以使用 workflow env create ``ENV_NAME`` 创建单个 conda 环境。

使用大量 conda 环境的 Snakemake 工作流程可能需要很长时间来安装，因为每个环境都是顺序创建的。运行 workflow env create --workers number_of_workers 将并行创建所有 conda 环境，最多可达请求的工作线程数（默认为 1）。

.. code-block:: bash

    chloroscan env create --workers 4  # create up to 4 conda envs at a time

.. warning::

    某些 conda 环境可能不支持并行创建。如果遇到错误，请尝试减少工作线程数。


env remove
-----------

此命令会删除工作流程中的所有 conda 环境。您也可以通过指定环境名称来删除单个环境。使用 ``--force`` 选项可以跳过确认提示。

.. code-block:: bash

    chloroscan env remove [OPTIONS] [ENV_NAME...]

``ENV_NAME``: 要删除的环境的名称。

env run
--------

env run 命令允许您在工作流程环境中运行命令。

.. code-block:: bash

    chloroscan env run --env ENV_NAME CMD...

``CMD...``: 要在指定环境中执行的命令及其参数。

确保将 ``ENV_NAME`` 替换为所需环境的实际名称，并将 ``CMD`` 替换为您要运行的命令。

例如：
.. code-block:: bash

    chloroscan env run -e my_environment "python script.py --input input_file.txt --output output_file.txt"

此命令将在工作流程中的 ``my_environment`` 环境中运行 ``python script.py --input input_file.txt --output output_file.txt`` 命令。
根据您的具体用例调整命令和环境名称。

env show
--------

此命令显示工作流程中使用的环境配置文件的内容。

.. code-block:: bash

    chloroscan env show [OPTIONS]


info
=========

info 子命令提供有关工作流程的 JSON 格式信息。

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

profile 子命令提供了多个命令来管理工作流程的配置文件。配置文件用于定义工作流程的不同配置，例如，您可以配置工作流程在 HPC 上的运行方式。您可以在 Snakemake 文档中了解有关配置文件的更多信息。

.. note

    为了使 snk 能够访问配置文件，配置文件必须位于工作流程的 profiles 目录中。

profile list
-------------

列出工作流程中的配置文件。
.. code-block:: bash

    chloroscan profile list

profile show
-------------

显示配置文件的内容。

.. code-block:: bash

    chloroscan profile show slurm

您可以通过将输出重定向到文件来保存配置文件。

.. code-block:: bash

    chloroscan profile show slurm > profile/config.yaml

.. note 

    您必须将配置文件保存为 config.yaml，以便工作流程可以访问它。

Load a profile
--------------

您可以使用 run 命令中的 ``--profile`` 选项加载配置文件。

.. code-block:: bash

    chloroscan run --profile profile/config.yaml

profile edit
--------------

edit 命令将在默认编辑器中打开配置文件。保存的更改将修改安装的默认配置文件设置。

.. code-block:: bash

    chloroscan profile edit slurm

github
==============

.. code-block:: bash

    chloroscan github

打开 ChloroScan GitHub 页面：https://github.com/andyargueasae/chloroscan

docs
==============

.. code-block:: bash

    chloroscan docs

打开 ChloroScan 文档页面：https://andyargueasae.github.io/chloroscan/
