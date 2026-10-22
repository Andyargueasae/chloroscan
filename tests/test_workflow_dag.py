"""Validate the full workflow DAG without downloading databases or running tools."""
from pathlib import Path
import subprocess


def test_full_workflow_dag(tmp_path):
    root = Path(__file__).resolve().parent
    data = root / "test-data"
    cmd = [
        "chloroscan", "run",
        "--Inputs-assembly", str(data / "test_plastid.fasta"),
        "--Inputs-depth-txt", str(data / "contig_depth_test.txt"),
        "--Inputs-batch-name", "CI_DRYRUN",
        "--outputdir", str(tmp_path / "output"),
        "--cat-database", str(data / "CAT_mini" / "db"),
        "--cat-taxonomy", str(data / "CAT_mini" / "tax"),
        "--dry", "-c", "1",
        "--directory", str(tmp_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=180)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert "This was a dry-run" in result.stdout + result.stderr
    assert "check_done" in result.stdout + result.stderr
