def test_feasibility(model):
    assert model.results.attrs["termination_condition"] == "optimal"
