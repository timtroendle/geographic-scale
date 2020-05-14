import pytest
import pandas as pd

DEMAND_TECH = "demand_elec"
EPSILON_IMPORTS = 0.001 # 0.1 %
EPSILON_BIOENERGY = 0.01 # 1 %


@pytest.fixture
def rel_net_import_threshold(model):
    autarky_scale, _, autarky_level, _, _ = model.inputs.scenario.split("-")
    share = 1 - float(autarky_level) / 100
    # if autarky_scale == "continental":
    #     share = 100 # basically unlimited national imports allowed
    return share


@pytest.fixture
def transmission_loc_techs_per_autarkic_group(model):
    return [
        getattr(model._model_data, f"group_constraint_loc_techs_{group}")
        for group in model._model_data.group_names_net_import_share_max.values
    ]


@pytest.fixture
def net_imports_per_autarkic_group(model, transmission_loc_techs_per_autarkic_group, scaling_factors):
    imports = [
        (model
         .results
         .carrier_prod
         .sel(loc_tech_carriers_prod=[f"{loc_tech.item()}::electricity" for loc_tech in transmission_loc_techs])
         .sum()
         .item()) / scaling_factors["power"]
        for transmission_loc_techs in transmission_loc_techs_per_autarkic_group
    ]
    exports = [
        (model
         .results
         .carrier_con
         .sel(loc_tech_carriers_con=[f"{loc_tech.item()}::electricity" for loc_tech in transmission_loc_techs])
         .sum()
         .item()) / scaling_factors["power"]
        for transmission_loc_techs in transmission_loc_techs_per_autarkic_group
    ]
    return [i + e for i, e in zip(imports, exports)]


@pytest.fixture
def demand_per_autarkic_group(carrier_con, transmission_loc_techs_per_autarkic_group):
    locs_per_autarkic_group = [
        list(set([loc_tech.item().split("::")[0] for loc_tech in transmission_loc_techs]))
        for transmission_loc_techs in transmission_loc_techs_per_autarkic_group
    ]
    return [
        carrier_con.sel(techs=DEMAND_TECH, locs=locs).sum().item()
        for locs in locs_per_autarkic_group
    ]


def test_biofuel_potential_not_exceeded(carrier_prod, biofuel_potential, location):
    biofuel_generation = carrier_prod.sel(techs="biofuel", locs=location).sum("timesteps").item()
    assert biofuel_generation <= biofuel_potential * (1 + EPSILON_BIOENERGY)


def test_net_imports(net_imports_per_autarkic_group, demand_per_autarkic_group, rel_net_import_threshold):
    rel_net_imports_per_autarkic_group = pd.Series([
        i / -d
        for i, d in zip(net_imports_per_autarkic_group, demand_per_autarkic_group)
    ])
    assert (rel_net_imports_per_autarkic_group <= (rel_net_import_threshold + EPSILON_IMPORTS)).all()
