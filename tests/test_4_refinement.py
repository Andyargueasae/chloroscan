def test_refinement(run_workflow):
    w = run_workflow("TEST_OUT/working/refined_bins")
    w.assert_exists()
    w.assert_exists("TEST_OUT/working/refinement_contig_summary.txt")