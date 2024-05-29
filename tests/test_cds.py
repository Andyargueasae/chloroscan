def test_cds(run_workflow):
    w = run_workflow("TEST_OUT/working/cds-extraction")
    w.assert_exists()
    w.assert_dir_exists("TEST_OUT/working/cds-extraction/GFFs")
    