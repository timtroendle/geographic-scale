import pytest

HYDRO_TECHS = ["pumped_hydro", "hydro_reservoir", "hydro_run_of_river"]


@pytest.mark.parametrize(
    "tech",
    HYDRO_TECHS
)
def test_no_hydro_costs(cost, tech):
    assert cost.sel(techs=tech).sum("locs") == pytest.approx(0.0)
