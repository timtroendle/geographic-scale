import pytest

HYDRO_TECHS = ["pumped_hydro", "hydro_reservoir", "hydro_run_of_river"]
CAPACITY_FACTORS = [# depends on capacity-factors/max and /min and time resolution configuration value
    ('hydro_run_of_river', 0.001 / 4, 1.0),
    ('hydro_reservoir', 0.001 / 4, 10.0),
    ('open_field_pv', 0.001 / 4, 1.0),
    ('roof_mounted_pv', 0.001 / 4, 1.0),
    ('wind_offshore', 0.001 / 4, 1.0),
    ('wind_onshore_competing', 0.001 / 4, 1.0),
    ('wind_onshore_monopoly', 0.001 / 4, 1.0),
]


@pytest.mark.parametrize(
    "tech",
    HYDRO_TECHS
)
def test_hydro_has_costs(cost, tech):
    assert cost.sel(techs=tech).sum("locs") > 0.0


@pytest.mark.parametrize(
    "tech,min_cf,max_cf",
    CAPACITY_FACTORS
)
def test_minimal_capacity_factors(resource, tech, min_cf, max_cf):
    cf = resource.sel(techs=tech)
    assert cf.min(["timesteps", "locs"]) >= 0
    assert cf.where(cf > 0).min(["timesteps", "locs"]) >= min_cf


@pytest.mark.parametrize(
    "tech,min_cf,max_cf",
    CAPACITY_FACTORS
)
def test_maximal_capacity_factor(resource, tech, min_cf, max_cf):
    assert resource.sel(techs=tech).min(["timesteps", "locs"]) <= max_cf
