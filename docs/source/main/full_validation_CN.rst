==========================
完整安装与运行验证
==========================

ChloroScan 的快速 CI 测试可以检查 Python 包安装、CLI 和 Snakemake 工作流 DAG。该测试**不会**下载大型数据库或执行完整的生物学分析。单独的可选测试可以创建工作流 Conda 环境，但这仍不代表整个分析已成功运行。

建议在 Linux x86_64、Python 3.10 环境中操作，准备足够的内存和磁盘空间。CAT 数据库解压后约需 85 GB。CORGI 默认使用 CPU；如需 CUDA 12.4 且设备兼容，请在创建 Conda 规则环境前设置 ``CHLOROSCAN_INSTALL_CUDA=1``。

1. 准备输入
------------

准备一个非空的 assembly FASTA 和对应 contig ID 的深度 TSV 文件，下载并解压 CAT 数据库，核实 ``db`` 目录存在且 ``tax/names.dmp`` 可读取。

2. 运行完整验证脚本
--------------------

克隆与当前发行版对应的 ChloroScan 仓库，在安装好 ChloroScan 和 Conda 的环境中，从仓库根目录运行：

.. code-block:: bash

    bash scripts/validate_full_install.sh /data/assembly.fasta \\
      /data/contig_depth.tsv sample_01 /data/chloroscan_results \\
      /data/CAT_database/db /data/CAT_database/tax 8

脚本将依次执行：检查工作流 DAG、创建所有 Snakemake Conda 规则环境、运行完整流程，并检查输出目录中的 ``ChloroScan.done``。

Binny 的 A2K 标记数据库将下载到 ``OUTPUT_DIR/resources/A2K_database``，无需写入已安装 Python 包的目录。

3. 常见故障排查
----------------

* **命令参数不识别**：运行 ``chloroscan run --help``；深度文件使用 ``--Inputs-depth-txt``，不是旧参数 ``--Inputs-depth-profile``。
* **CAT 路径错误**：确认 ``--cat-database`` 指向数据库目录，``--cat-taxonomy`` 指向包含 ``names.dmp`` 的目录。
* **Conda 环境创建失败**：使用 ``--use-conda --conda-frontend conda``，检查网络连接、磁盘空间、Conda 缓存和软件源；不必强制降级旧版 Conda/Mamba。
* **CUDA/PyTorch 不兼容**：先使用默认 CPU 环境；仅在兼容 GPU 和驱动上启用 ``CHLOROSCAN_INSTALL_CUDA=1``。
* **工作流中途失败**：检查输出目录 ``logging_info/`` 内的日志以及 Snakemake 提示，确认输入 contig ID 与深度表一致、数据库可读且输出目录可写。

GitHub Actions 中 ``testing`` 工作流的可选 ``provision_envs`` 检查仅验证环境创建。要确认实际分析可用，请在具备数据库及资源的工作站或 HPC 节点上完成全流程。