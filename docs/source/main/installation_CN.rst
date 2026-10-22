====================
安装 ChloroScan
====================

ChloroScan 持续更新。建议在 Linux x86_64 系统上使用 Python 3.10，并先通过小型测试数据验证安装。Python 包元数据允许 Python 3.9 及以上版本，但完整工作流尚未覆盖所有 Python 版本和操作系统。

1. 安装 Python 和 Conda
-----------------------

推荐安装 Miniforge：https://github.com/conda-forge/miniforge/releases/ 。下载适合您系统的安装文件（例如 ``Miniforge3-Linux-x86_64.sh``），然后执行：

.. code-block:: bash

    bash Miniforge3-Linux-x86_64.sh
    conda init bash
    source ~/.bashrc
    conda --version

对于 Snakemake 7.32，指定 ``--conda-frontend conda``；无需降级到 Conda 23.3.1 或 Mamba 1.4.2。

创建独立的 Python 环境并安装最新发行版：

.. code-block:: bash

    conda create -n chloroscan_env python=3.10 -y
    conda activate chloroscan_env
    python -m pip install chloroscan
    chloroscan --help
    chloroscan run --help

亦可使用 ``python3.10 -m venv chloroscan_env`` 和 ``source chloroscan_env/bin/activate`` 创建虚拟环境，前提是 Conda 仍可用于 Snakemake 每个规则的环境。

对于编码序列提取，项目提供了 FragGeneScanRs 安装辅助脚本，请从项目仓库下载并检查其先决条件：``install_fraggenescanrs.sh``。

2. CAT 数据库
--------------

运行时必须明确提供 CAT 的数据库目录 ``--cat-database`` 和分类目录 ``--cat-taxonomy``。推荐使用项目整理的质体蛋白数据库：https://doi.org/10.26188/27990278 。

.. code-block:: bash

    figshare download -o CAT_db.tar.gz 27990278
    mkdir -p /path/to/CAT_database
    tar -xzf CAT_db.tar.gz -C /path/to/CAT_database

解压后请检查数据库文件的实际目录层级；``tax`` 目录应包含 ``names.dmp``。数据库压缩包约 47 GB，解压后约 85 GB，应预留足够磁盘空间。

3. 配置 Snakemake 环境
-----------------------

如已有排序的 BAM 比对文件目录，可使用 ``--Inputs-alignment``。如果已经生成了 contig 测序深度 TSV 文件，则改用 ``--Inputs-depth-txt``；无需同时提供两者。

以下示例只创建规则环境，不运行完整分析：

.. code-block:: bash

    chloroscan run --Inputs-assembly /data/assembly.fasta \
      --Inputs-depth-txt /data/contig_depth.tsv \
      --Inputs-batch-name sample_01 --outputdir /data/chloroscan_output \
      --cat-database /data/CAT_database/db --cat-taxonomy /data/CAT_database/tax \
      --use-conda --conda-frontend conda --conda-create-envs-only \
      --conda-prefix /data/chloroscan_conda_envs --cores 8

运行 ChloroScan 时，Binny 标记数据库会下载到 ``OUTPUT_DIR/resources/A2K_database``，不是在创建 Conda 环境时下载。请确保输出目录可写且空间充足。

4. CPU/GPU 选择
----------------

CORGI 默认使用 CPU 兼容环境。仅在具有兼容 NVIDIA 驱动和 GPU 的设备上，如确实需要 CUDA 12.4 加速，可在创建规则环境**之前**设置：

.. code-block:: bash

    export CHLOROSCAN_INSTALL_CUDA=1

不具备兼容 GPU 时不要设置该变量。

5. 完整测试和排错
-----------------

使用 ``chloroscan run --help`` 检查当前安装版本的真实选项。请确认输入文件非空、CAT 路径有效、``tax/names.dmp`` 存在且磁盘空间足够。完整执行后检查 ``OUTPUT_DIR/ChloroScan.done``。

请阅读 :doc:`full_validation_CN`，了解首次运行的测试脚本、数据准备和常见故障。