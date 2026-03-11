============
下载
============

我们仍在升级 ChloroScan 以更好地还原环境中的叶绿体基因组MAG，为了获得最佳效果，请确保安装最新版本。

1. pip 和 conda/mamba
----------------------

要安装工作流程，请使用 pip。首先创建一个虚拟环境。确保您的计算机严格遵循 Python 版本限制：> 3.9 且 < 3.12。作为参考，我们用Python 3.10.4 开发了ChloroScan。

运行此命令以检查您的 Python 版本是否兼容：

.. code-block:: bash

    python3 --version

另一个关键的包是mamba，它是您能够安装ChloroScan工作流中每一个步骤的虚拟环境的关键。您可以通过从以下链接下载 miniforge 来安装 mamba：

.. code-block:: bash

    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh


下载完成后，运行以下命令安装 mamba（将 {OPERATING_SYSTEM} 和 {ARCHITECTURE} 替换为您自己的系统信息，例如 Miniforge3-Linux-x86_64.sh）：

.. code-block:: bash

    bash Miniforge3-{OPERATING_SYSTEM}-{ARCHITECTURE}.sh

可以通过运行以下命令检查 conda 和 mamba 是否完整安装：

.. code-block:: bash

    # 为 bash 初始化 conda 和 mamba。
    conda init bash
    mamba init bash
    source ~/.bashrc

    conda --version
    mamba --version

    # 查询 mamba 和 conda 的安装路径
    which mamba
    which conda

``注意``: Mamba 和 conda 的版本会直接影响安装，作为参考，我们在开发ChloroScan的时候用的是conda=23.3.1和mamba=1.4.2。

运行以下命令以设置虚拟环境并安装 ChloroScan。

.. code-block:: bash

    python3 -m venv chloroscan_env
    source chloroscan_env/bin/activate

最后确定以上条件达成后，运行以下命令安装 ChloroScan：

.. code-block:: bash

    pip3 install chloroscan

或者，您可以创建一个具有兼容 Python 版本的 conda 环境，并在其中安装 chloroscan：

.. code-block:: bash

    mamba env create -n chloroscan_env python=3.10
    conda activate chloroscan_env
    which pip; which python; # 检查 pip 和 python 是否来自 conda 环境，这是必须保证的。 
    #!注意!: 系统的默认 Python （/usr/bin/python） 不应该显示在这里，否则将会导致令人头疼的问题。
    pip3 install chloroscan
    pip3 install chloroscan
    
除了安装，ChloroScan 还需要几个数据库：一个用于 CAT 的分类学分类，另一个用于通过分箱推断质体 MAG 的质量。

2. CAT 数据库
-------------------

您的选择可以是使用包含 NCBI nr 数据库中所有蛋白质的 CAT 预处理数据库。然而，在运行 CAT 时，它占用大量内存且耗时，并且需要约 300GB 的磁盘空间。

为了评估叶绿体 contig 的物种分类，我们建议下载我们精心策划的蛋白质数据库，该数据库包含来自 NCBI 的所有叶绿体基因组编码蛋白质和来自 UniRef90 的蛋白质（解压后大小约为 85GB）。

下载命令：

.. code-block:: bash

   figshare download -o CAT_db.tar.gz 27990278

3. Binny 数据库和 conda 环境
------------------------------

您无需手动下载 binny 的数据库，ChloroScan 会在为每个作业安装 conda 环境后自动加载它。

运行 chloroscan 并添加标志 "``--conda-create-envs-only``" 以设置 conda 环境。记得始终添加 "``--use-conda``" 标志以确保工作流程使用 conda 环境。您可以使用 "``--conda-prefix``" 标志指定 conda 环境的前缀，否则环境将创建在默认位置。

先版本有四个参数需要一直添加： ``Inputs-assembly``、 ``Inputs-batch-name``、 ``outputdir`` 和 ``Inputs-alignment``。（您可以使用测试数据（test data）来先进行conda环境配置）
.. code-block:: bash

    chloroscan run --use-conda --conda-create-envs-only --conda-prefix="/path/to/your/conda/envs" \
        --Inputs-assembly "/path/to/your/assembly.fasta" --Inputs-batch-name "sample_name" \
        --outputdir "/path/to/your/output_dir" --Inputs-alignment "/path/to/your/bams_directory" \

坐到这儿就标志着 ChloroScan 的准备工作已经完成了，您可以开始运行它。 
