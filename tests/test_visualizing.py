import os
def test_visualizing(run_workflow):
    w = run_workflow("TEST_OUT/working/visualizations")
    w.assert_exists()
    w.assert_exists("TEST_OUT/working/visualizations/Scatter_GCLogDepth.png")
    w.assert_exists("TEST_OUT/working/visualizations/LogDepth_Violin.png")
    assert len([figure for figure in os.listdir("TEST_OUT/working/visualizations") if figure.endswith("_taxonomy_composition.png")]) == 1, "There should be exactly one visualization output directory"
