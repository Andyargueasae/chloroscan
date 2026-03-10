
======================
运行 ChloroScan
======================

ChloroScan 专注于叶绿体基因组 MAGs 的恢复，特别是 ``光合微藻和大型藻类，及其衍生的复杂质体，例如硅藻和定鞭藻``。它主要需要的输入是 contigs 及其测序深度。
在大多数情况下，默认设置就能很好地工作，但你也可以为每个步骤配置参数。在这里，我们将简要介绍每个步骤的输入和输出，以及你可以为每个步骤配置的参数。

1. 数据中是否存在叶绿体？
======================================

光合藻类/光合原生生物通常栖息在海洋水域、淡水、地衣和潮湿土壤等环境中。人为环境、地下环境和沉积物中它们的丰度通常较低。因此，如果你的数据来自这些环境，ChloroScan 可能无法恢复任何叶绿体 MAGs。

同时，如果你的数据包含过于碎片化的叶绿体 contigs 且覆盖度较低（通常 < 5.0 x），ChloroScan 仍可能无法检测到它们。希望我们更新的版本能够解决这些问题。

2. 常用的最简命令
==========================

你已经有了 contigs、bam 格式的映射文件和配置文件。现在你可以用一个命令运行整个工作流程：

.. code-block:: bash

   chloroscan run --Inputs-assembly input_contigs.fasta --Inputs-alignment PATH/to/bams \
      --Inputs-batch-name "my_batch" --outputdir Path/to/output --use-conda --cores=12 \
      --cat-database PATH/to/CAT_db/db --cat-taxonomy path/to/CAT_db/tax

有时你可能还希望获得以表格形式表示的深度信息。它看起来像这样：

.. code-block::
   
   S0C861	3.52491757
   S0C1664	2.73124830
   S0C2713	12.64139886
   S0C3242	2.51473363
   S0C8106	23.82718202
   S0C8631	2.69335600
   S0C9609	2.49439900
   ...

第一列是 contig 的 ID，第二列是 contig 的平均深度。我们也接受这种格式，只需将命令更改为：

.. code-block:: bash

   DEPTH_PROFILE=path/to/depth_profile.tsv
   chloroscan run --Inputs-assembly input_contigs.fasta --Inputs-depth-profile $DEPTH_PROFILE \
      --Inputs-batch-name "my_batch" --outputdir Path/to/output --use-conda --cores=12 \
      --cat-database PATH/to/CAT_db/db --cat-taxonomy path/to/CAT_db/tax

有关输入的更多详细信息，请查看 ``输入和输出`` 部分。
如果你想使用我们的测试数据，可以使用 ``README`` 中显示的命令下载测试数据：
  
.. code-block:: bash
   
   figshare download -o simulated_metagenomes.tar.gz 28748540


3. 命令参数说明
========================================

ChloroScan 的完整命令空间如下所示：

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

以下列出了 ChloroScan 的所有参数。
 - ``--Inputs-assembly``: Fasta 格式的 contigs 装配文件路径，包含各种生物的 contigs。
 - ``--Inputs-depth-txt``: 样本中每个 contig 丰度的制表符分隔文本文件路径。第一列是 contig 的 ID，第二列是 contig 的平均深度。我们也接受这种格式，只需将命令更改为 ``--Inputs-depth-profile``。
 - ``--Inputs-alignment``: 包含 contigs 对齐文件的文件夹路径。对齐文件应为 bam 格式，并命名为 "sample_name.bam"。样本名称将从 bam 文件名中提取，通过去掉 ".bam" 后缀。样本名称将在下游分析和输出文件中使用。
 - ``--Inputs-batch-name``: 批次名称。该名称将在下游分析和输出文件中使用。
 - ``--outputdir``: 工作流的输出目录路径。最终结果将存储在此目录中。中间结果将存储在输出目录下名为 "working" 的子目录中。默认值为 "output"。
 - ``--tmpdir``: 工作流的临时目录路径。默认值为 "tmp"。
 - ``--binning-universal-length-cutoff``: 通用分箱的长度阈值。短于此长度的 contigs 将在分箱前被过滤。默认值为 1500bp。
 - ``--binning-snakemake-env``: 用于运行 binny 的自定义 snakemake 环境。如果未指定，将使用默认的 conda 环境。
 - ``--binning-mantis-env``: 用于安装 mantis_pfa 并注释基因的自定义 Mantis 虚拟环境。如果未指定，将使用默认的 conda 环境。
 - ``--binning-outputdir``: 分箱的输出目录路径。默认值为 "binny_output"。
 - ``--binning-clustering-epsilon-range``: HDBSCAN 聚类中重要的 epsilon的取值范围。默认值为 "0.250,0.000"。
 - ``--binning-clustering-hdbscan-min-sample-range``: HDBSCAN 聚类中 min_samples 的取值范围，较大的值意味着较大的 MAGs。默认值为 "1,4,7,10"。
 - ``--binning-bin-quality-purity``: 分箱质量的最小纯度。默认值为 95。
 - ``--binning-bin-quality-starting-completeness``: 分箱质量的起始完整性。Binny 使用滑动完整性来过滤分箱。默认值为 92.5。
 - ``--binning-bin-quality-min-completeness``: 分箱质量的最小完整性。默认值为 72.5。
 - ``--corgi-min-length``: CORGI 处理的 contigs 的最小长度。默认值为 500bp。
 - ``--corgi-save-filter``: 保存 CORGI 过滤后的 contigs（注意：可能需要较长时间）。默认值为 no-corgi-save-filter。
 - ``--corgi-batch-size``: CORGI 推断 contigs 分类标签的批次大小。默认值为 1。
 - ``--corgi-pthreshold``: CORGI 判断 contigs 类别是否为真正的叶绿体或其他类别的 P 值阈值。默认值为 0.9。
 - ``--cat-database``: CAT 数据库的路径，保存经过 diamond 处理的蛋白质序列。默认值为 "PATH/TO/CAT_db/db"。
 - ``--cat-taxonomy``: CAT 数据库的分类标签路径。默认值为 "PATH/TO/CAT_db/tax"。
 - ``--krona-env``: Krona 环境的路径。默认值为 "kronatools"。 **现在通常不需要设置此项**。