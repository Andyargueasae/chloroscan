
# Only you might fail, what is the exact reason?
def test_binny(run_workflow):
    w = run_workflow("TEST_OUT/working/binny")
    w.assert_exists()
