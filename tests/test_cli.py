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
    # breakpoint()
    print(result.stdout)
    print(result.stderr)
    assert result.exit_code == 0
    assert result.stderr == ""
    assert "Path to fasta format assembly of contigs from all sorts of organisms." in result.stdout
    assert "Path to a tab-separated text storing abundance of each contig in the sample." in result.stdout
    assert "Path to the folder containing alignment files of the contigs." in result.stdout
    
        