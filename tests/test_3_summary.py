
# Only you might fail, what is the exact reason?
def test_summary(run_workflow):
    w = run_workflow("TEST_OUT/working/summary/cross_ref.tsv")
    w.assert_exists()
    w.assert_contains("k127_5658")
    w.assert_contains("k127_37878")
    w.assert_contains("k127_14538")
    w.assert_contains("k127_42002")

