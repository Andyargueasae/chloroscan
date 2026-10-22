#!/usr/bin/env bash
# Prepare, validate, and run ChloroScan with user-owned resource directories.
set -euo pipefail

if [[ "$#" -lt 6 || "$#" -gt 8 ]]; then
  echo "Usage: $0 ASSEMBLY_FASTA DEPTH_TSV BATCH_NAME OUTPUT_DIR CAT_DB CAT_TAX [THREADS=8] [CONDA_PREFIX=OUTPUT_DIR/conda_envs]" >&2
  exit 2
fi

assembly=$1
depth=$2
batch=$3
out=$4
cat_db=$5
cat_tax=$6
threads=${7:-8}
conda_prefix=${8:-$out/conda_envs}

command -v chloroscan >/dev/null || { echo "chloroscan is missing from PATH" >&2; exit 1; }
command -v conda >/dev/null || { echo "conda is missing from PATH" >&2; exit 1; }
[[ -s "$assembly" ]] || { echo "Assembly is missing or empty: $assembly" >&2; exit 1; }
[[ -s "$depth" ]] || { echo "Depth table is missing or empty: $depth" >&2; exit 1; }
[[ -d "$cat_db" ]] || { echo "CAT database directory is missing: $cat_db" >&2; exit 1; }
[[ -f "$cat_tax/names.dmp" ]] || { echo "CAT taxonomy names.dmp is missing: $cat_tax" >&2; exit 1; }

mkdir -p "$out" "$conda_prefix"
out=$(cd "$out" && pwd)
conda_prefix=$(cd "$conda_prefix" && pwd)

args=(
  --Inputs-assembly "$assembly"
  --Inputs-depth-txt "$depth"
  --Inputs-batch-name "$batch"
  --outputdir "$out"
  --cat-database "$cat_db"
  --cat-taxonomy "$cat_tax"
  --use-conda
  --conda-frontend=conda
  --conda-prefix="$conda_prefix"
  -c "$threads"
)

echo "[1/3] Checking the complete workflow DAG"
chloroscan run "${args[@]}" --dry

echo "[2/3] Creating all Conda rule environments"
chloroscan run "${args[@]}" --conda-create-envs-only

echo "[3/3] Running the complete workflow"
chloroscan run "${args[@]}"

test -f "$out/ChloroScan.done"
echo "ChloroScan completed successfully: $out/ChloroScan.done"
