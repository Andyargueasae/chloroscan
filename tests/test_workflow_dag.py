"""Validate workflow DAGs without downloading databases or running tools."""
import re
import subprocess
from pathlib import Path

import pytest


DATA = Path(__file__).resolve().parent / "test-data"
CAT_ARGS = [
    "--cat-database", str(DATA / "CAT_mini" / "db"),
    "--cat-taxonomy", str(DATA / "CAT_mini" / "tax"),
]
CAT_RULES = {
    "CAT_taxonomy_identification_plus_annotation",
    "documenting_binny_results",
    "refine_bins",
    "visualize_results",
    "krona_taxonomy_plot",
}


def dry_run(tmp_path, *args):
    cmd = [
        "chloroscan", "run",
        "--Inputs-assembly", str(DATA / "test_plastid.fasta"),
        "--Inputs-batch-name", "CI_DRYRUN",
        "--outputdir", str(tmp_path / "output"),
        "--dry", "-c", "1",
        "--directory", str(tmp_path),
        *args,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=180)
    return result.returncode, result.stdout + "\n" + result.stderr


@pytest.mark.parametrize("mode", ["default", "cat_paths_only", "use_cat"])
@pytest.mark.parametrize("coverage", ["depth_table", "alignment"])
def test_workflow_dag(tmp_path, mode, coverage):
    args = []
    if coverage == "alignment":
        bams = tmp_path / "bams"
        bams.mkdir()
        args.extend(["--Inputs-alignment", str(bams)])
    else:
        args.extend(["--Inputs-depth-txt", str(DATA / "contig_depth_test.txt")])

    if mode != "default":
        args.extend(CAT_ARGS)
    if mode == "use_cat":
        args.append("--use-cat")

    code, output = dry_run(tmp_path, *args)
    assert code == 0, output
    assert "This was a dry-run" in output
    job_stats = output.split("Job stats:", 1)[1].split("\n\n", 1)[0]
    rules = set(re.findall(r"^(\w+)[ \t]+[1-9]\d*$", job_stats, re.MULTILINE))
    assert {"corgi_prediction", "binny_run", "check_binny_output", "cds_extraction", "check_done"} <= rules
    expected_coverage = {"contig_depth_calculation", "merge_contig_depth"} if coverage == "alignment" else set()
    assert {"contig_depth_calculation", "merge_contig_depth"} & rules == expected_coverage
    assert CAT_RULES & rules == (CAT_RULES if mode == "use_cat" else set())


@pytest.mark.parametrize("cat_args", [CAT_ARGS[:2], CAT_ARGS[2:]])
def test_use_cat_requires_both_paths(tmp_path, cat_args):
    code, output = dry_run(
        tmp_path, "--Inputs-depth-txt", str(DATA / "contig_depth_test.txt"),
        "--use-cat", *cat_args,
    )
    assert code != 0, output
    assert "--use-cat requires both --cat-database and --cat-taxonomy" in output
