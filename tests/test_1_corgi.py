import pytest

@pytest.mark.dependency()
def test_corgi(run_workflow):
    w = run_workflow("TEST_OUT/working/corgi/plastid.fasta")
    w.assert_exists()
    w.assert_contains(">k127_14538")
    w.assert_contains(">k127_37878")
    w.assert_contains(">k127_5658")
    w.assert_contains(">k127_42002")