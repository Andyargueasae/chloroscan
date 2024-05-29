def test_visualizing(run_workflow):
    w = run_workflow("TEST_OUT/working/visualizations")
    w.assert_exists()