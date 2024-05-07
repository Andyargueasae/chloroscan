# Next step: test those rules.
from snk_cli.testing import SnkCliRunner
from chloroscan.cli import chloroscan



def test_help():
    runner = SnkCliRunner(chloroscan)
    result = runner.invoke(["--help"])
    assert result.exit_code == 0
    assert result.stderr == ""
    assert "A state-of-art metagenomic workflow extracting chloroplast genomes" in result.stdout
    

def test_run_help():
    runner = SnkCliRunner(chloroscan)
    result = runner.invoke(["run", "--help"])
    assert result.exit_code == 0
    assert result.stderr == ""
    assert "Path to fasta format assembly of contigs " in result.stdout
        