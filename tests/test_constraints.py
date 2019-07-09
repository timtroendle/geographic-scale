import pytest

DEMAND_TECH = "demand_elec"


@pytest.fixture
def net_import_threshold(carrier_con, model):
    demand = carrier_con.sel(techs=DEMAND_TECH).sum(["locs", "timesteps"]).item()
    autarky_scale, _, autarky_level, _, _ = model.inputs.scenario.split("-")
    share = 1 - float(autarky_level) / 100
    if autarky_scale == "continental":
        share = 100 # basically unlimited national imports allowed
    return demand * -1 * share


@pytest.fixture
def net_imports(model, transmission_carrier_prod, units):
    if "loc_techs_transmission" in model.inputs:
        edges = [(line.split(":")[0], line.split(":")[-1]) for line in model.inputs.loc_techs_transmission.values]
    else:
        edges = []
    edges = filter( # transform to links between countries
        lambda x: x[0] != x[1],
        map(
            lambda x: (units.loc[x[0], "country_code"], units.loc[x[1], "country_code"]),
            edges
        )
    )
    edges = list(set(tuple(sorted((a, b))) for a, b in edges)) # filter duplicates
    return sum([
        abs(link_timeseries(transmission_carrier_prod, units, source, sink).sum("timesteps"))
        for (source, sink) in edges
    ])


def link_timeseries(transmission_carrier_prod, units, country_code_a, country_code_b):
    gen = transmission_carrier_prod
    nodes_a = units[units.country_code == country_code_a].index
    nodes_b = units[units.country_code == country_code_b].index
    a_to_b = gen.sel(source=nodes_a, sink=nodes_b).sum(["sink", "source"])
    b_to_a = gen.sel(source=nodes_b, sink=nodes_a).sum(["sink", "source"])
    return a_to_b - b_to_a


@pytest.fixture
def transmission_carrier_prod(carrier_prod):
    gen = carrier_prod.sel(techs=carrier_prod.techs.str.contains("ac_transmission"))
    return (gen.assign_coords(techs=[tech.item().split(":")[-1] for tech in gen.techs])
               .rename(locs="sink", techs="source"))


def test_biofuel_potential_not_exceeded(carrier_prod, biofuel_potential, location):
    biofuel_generation = carrier_prod.sel(techs="biofuel", locs=location).sum("timesteps").item()
    assert biofuel_generation <= biofuel_potential


@pytest.mark.xfail(reason="losses can currently be fed from net imports, see Calliope #270")
def test_net_imports(net_imports, net_import_threshold):
    # this tests tests imports on the national level only
    # thereby continental autarky will always pass which is fine
    # regional autarky will pass as long as the weaker national autarky is hold
    assert net_imports <= net_import_threshold
