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
    assert "--Inputs-assembly" in result.stdout
    assert "--Inputs-depth-txt" in result.stdout
    assert "--Inputs-alignment" in result.stdout
    assert "--Inputs-batch-name" in result.stdout
    assert "--outputdir" in result.stdout
    assert "--cat-database" in result.stdout
    assert "--cat-taxonomy" in result.stdout
    assert "--krona-env" in result.stdout
    assert "--binning" in result.stdout
    assert "--corgi" in result.stdout
    assert "--tmpdir" in result.stdout
    
        