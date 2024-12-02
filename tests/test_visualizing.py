def test_visualizing(run_workflow):
    w = run_workflow("TEST_OUT/working/visualizations")
    w.assert_exists()
    w.assert_exists("TEST_OUT/working/visualizations/Scatter_GCLogDepth.png")
    w.assert_exists("TEST_OUT/working/visualizations/LogDepth_Violin.png")
    w.assert_exists("TEST_OUT/working/visualizations/TEST_OUT_binny_I04R01.000000_C84_P97_Haptophyta_taxonomy_composition.png")
