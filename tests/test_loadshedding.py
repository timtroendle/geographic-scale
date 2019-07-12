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


@pytest.mark.xfail(reason="Load shedding can and is currently exported.") # FIXME don't export load shedding
def test_loadshedding_not_exported(carrier_prod, carrier_con, location_with_loadshedding):
    load_shedding = carrier_prod.sel(locs=location_with_loadshedding, techs=LOAD_SHEDDING)
    demand = carrier_con.sel(locs=location_with_loadshedding, techs=DEMAND)
    assert ((load_shedding + demand) <= 0).all()
