def test_krona(run_workflow):
    w = run_workflow("TEST_OUT/Krona.html")
    w.assert_exists()