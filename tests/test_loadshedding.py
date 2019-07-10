import pytest

LOAD_SHEDDING = "load_shedding"
DEMAND = "demand_elec"


@pytest.fixture
def location_with_loadshedding(carrier_prod, location):
    if LOAD_SHEDDING not in carrier_prod.techs:
        pytest.skip(f"Location {location} has no load shedding.")
    else:
        load_shedding = carrier_prod.sel(locs=location, techs=LOAD_SHEDDING)
        if not load_shedding.to_series().any():
            pytest.skip(f"Location {location} has no load shedding.")
        else:
            return location


def test_loadshedding_not_exported(carrier_prod, carrier_con, location_with_loadshedding):
    load_shedding = carrier_prod.sel(locs=location_with_loadshedding, techs=LOAD_SHEDDING)
    max_load_shedding_power = load_shedding.max("timesteps")
    peak_demand = carrier_con.sel(locs=location_with_loadshedding, techs=DEMAND).min("timesteps")
    assert peak_demand < 0
    assert max_load_shedding_power < peak_demand * -1
