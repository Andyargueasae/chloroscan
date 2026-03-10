============
下载
============

我们仍在升级 ChloroScan 以更好地还原环境中的叶绿体基因组MAG，为了获得最佳效果，请确保安装最新版本。

要安装工作流程，请使用 pip。首先创建一个虚拟环境。确保您的计算机严格遵循 Python 版本限制：> 3.9。

运行此命令以检查您的 Python 版本是否兼容：
.. code-block:: bash

    python3 --version

如果您的 Python 版本在 3.9 和 4.0 之间，您将能够继续设置虚拟环境并安装 ChloroScan。
运行以下命令以设置虚拟环境并安装 ChloroScan。

.. code-block:: bash

    python3 -m venv chloroscan_env
    source chloroscan_env/bin/activate

最后确定以上条件达成后，运行以下命令安装 ChloroScan：
.. code-block:: bash

    pip3 install chloroscan

注意，我们没有将 ChloroScan 上传到 conda。

除了安装，ChloroScan 还需要几个数据库：一个用于 CAT 的分类学分类，另一个用于通过分箱推断质体 MAG 的质量。

CAT 数据库
-------------------

您的选择可以是使用包含 NCBI nr 数据库中所有蛋白质的 CAT 预处理数据库。然而，在运行 CAT 时，它占用大量内存且耗时，并且需要约 300GB 的磁盘空间。

为了评估叶绿体 contig 的物种分类，我们建议下载我们精心策划的蛋白质数据库，该数据库包含来自 NCBI 的所有叶绿体基因组编码蛋白质和来自 UniRef90 的蛋白质（解压后大小约为 85GB）。

下载命令：

.. code-block:: bash

   figshare download -o CAT_db.tar.gz 27990278

Binny 数据库和 conda 环境
-------------------

您无需手动下载 binny 的数据库，ChloroScan 会在为每个作业安装 conda 环境后自动加载它。

运行 chloroscan 并添加标志 "``--conda-create-envs-only``" 以设置 conda 环境。
.. code-block:: bash

    chloroscan run --conda-create-envs-only

坐到这儿就标志着 ChloroScan 的准备工作已经完成了，您可以开始运行它。 
