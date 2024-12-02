def test_binny(run_workflow):
    w = run_workflow("TEST_OUT/working/CAT")
    w.assert_exists()
    w.assert_dir_exists("TEST_OUT/working/CAT")
    w.assert_exists("TEST_OUT/working/CAT/out.CAT.contig2classification.txt")
    # w.assert_include("bins/*.fasta")
    # Perhaps add new outputs to the rule and test whether there are. 
