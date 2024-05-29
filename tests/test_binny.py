
# Only you might fail, what is the exact reason?
def test_binny(run_workflow):
    w = run_workflow("TEST_OUT/working/binny")
    w.assert_exists()
    w.assert_dir_exists("TEST_OUT/working/binny/bins")
    # w.assert_include("bins/*.fasta")
    # Perhaps add new outputs to the rule and test whether there are. 

