

# ChloroScan

**A computational workflow designed to recover plastid genomes from metagenomes.**

ChloroScan is a Snakemake-based bioinformatics pipeline that integrates advanced filtering, binning, and taxonomic profiling tools to efficiently extract and assemble chloroplast (plastid) genomes from complex metagenomic datasets.

## 🌟 Features
- **Automated Plastid Identification**: Utilizes **CORGI** to rapidly filter and classify chloroplast contigs from whole metagenomic assemblies.
- **Robust Binning Workflow**: Leverages a customized **Binny** pipeline with HDBSCAN clustering to group contigs into high-quality Metagenome-Assembled Organisms (MAGs).
- **Taxonomic Profiling**: Integrates **CAT** (Composition-based Assignment Tool) for accurate, homology-based taxonomy prediction.
- **Highly Configurable**: Tunable parameters for contig length cutoffs, clustering epsilon ranges, completeness/purity thresholds, and batch processing.
- **Reproducible & Scalable**: Built on Snakemake, optimized for local machines and High-Performance Computing (HPC) clusters.

## 📦 Installation

ChloroScan requires Python `>=3.9`. It is recommended to install using `Poetry` or `pip`.

```bash
# Clone the repository
git clone https://github.com/Andyargueasae/chloroscan.git
cd chloroscan

# Install dependencies (Poetry recommended)
pip install poetry
poetry install

# Alternatively, install directly via pip
pip install .
```

> **💡 Note**: The binning module relies on Snakemake and Conda/Mamba environments. Ensure you have `snakemake` and `mamba` installed in your environment for full workflow functionality.

## 🚀 Usage

Execute the pipeline via the command line using the `chloroscan run` command.

### Basic Example
If you have a contig depth file available:
```bash
chloroscan run \
  --Inputs-assembly-path assembly.fasta \
  --Inputs-depth-txt contig_depth.txt \
  --Inputs-batch-name SAMPLE_01 \
  --outputdir chloroscan_results \
  --cores 12
```

### Without Depth File
If depth data is unavailable, you can provide the alignment folder containing sorted `.bam` files:
```bash
chloroscan run \
  --Inputs-assembly-path assembly.fasta \
  --Inputs-batch-name SAMPLE_01 \
  --outputdir chloroscan_results \
  --alignment-folder path/to/alignment_bams \
  --cores 12
```

### 🔧 Command Line Arguments
| Flag | Description | Required | Default |
|---|---|---|---|
| `--Inputs-assembly-path` | Path to FASTA assembly of contigs | ✅ | - |
| `--Inputs-depth-txt` | Path to tab-separated contig abundance/depth file | ❌ | - |
| `--alignment-folder` | Path to folder containing `.bam` read alignments | ❌ | - |
| `--Inputs-batch-name` | Unique identifier for the sample/batch | ✅ | - |
| `--outputdir` | Directory for final workflow outputs | ❌ | `TEST_OUT` |
| `--tmpdir` | Directory for intermediate temporary files | ❌ | `tmp` |
| `--cores` | Number of CPU threads to utilize | ❌ | System default |

### ⚙️ Advanced Configuration
You can fine-tune pipeline behavior by adjusting CLI flags or modifying the YAML configuration files in `chloroscan/config/`:

- **CORGI Filtering**: `--corgi-settings-min-length` (default: 1000), `--corgi-settings-pthreshold` (default: 0.50), `--corgi-settings-batch-size`
- **Binning Parameters**: `--binning-universal-length-cutoff`, `--binning-clustering-epsilon_range`, `--binning-bin-quality-min_completeness`
- **CAT Databases**: `--cat-database`, `--cat-taxonomy` (Ensure absolute paths are used for external databases)

## 📖 Documentation
For comprehensive guides, a beginner's tutorial, and detailed parameter explanations, visit the official documentation:
👉 [ChloroScan Documentation](https://andyargueasae.github.io/chloroscan/)

## 👥 Credits
Developed by:
- Andy Tong
- Robert Turnbull
- Vanessa Rossetto Marcelino
- Heroen Verbruggen

## 📜 License
This project is licensed under the **Apache-2.0 License**. See the `LICENSE` file for details.
