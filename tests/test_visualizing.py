def test_visualizing(run_workflow):
    w = run_workflow("TEST_OUT/working/visualizations")
    w.assert_exists()
    w.assert_exists("TEST_OUT/working/visualizations/Scatter_GCLogDepth.png")
    w.assert_exists("TEST_OUT/working/visualizations/LogDepth_Violin.png")
    w.assert_exists("TEST_OUT/working/visualizations/*_Haptophyta_taxonomy_composition.png")