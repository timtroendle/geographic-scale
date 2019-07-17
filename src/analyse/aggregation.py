"""Collect globally aggregated results."""
import functools
from pathlib import Path
from dataclasses import dataclass

import calliope
import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx


CARRIER = "electricity"
ONSHORE_WIND_TECHS = ["wind_onshore_competing", "wind_onshore_monopoly"]
PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
RE_TECHS = ["open_field_pv", "roof_mounted_pv", "wind_onshore_monopoly",
            "wind_onshore_competing", "wind_offshore", "hydro_run_of_river"]
ELECTRICITY_DEMAND_TECH = "demand_elec"
CACHE_SIZE = 20
M_TO_KM = 1e-3


@dataclass
class Variable:
    name: str
    unit: str
    value_function: callable = None


def excavate_all_results(paths_to_scenarios, path_to_units, scaling_factors, path_to_output):
    """Collect globally aggregated results."""
    scenarios = {
        Path(path_to_scenario).parent.name: calliope.read_netcdf(path_to_scenario)
        for path_to_scenario in paths_to_scenarios
    }
    units = (gpd.read_file(path_to_units)
                .set_index("id")
                .rename(index=lambda idx: idx.replace(".", "-"))
                .rename_axis(index="locs"))
    variables = _set_up_variables(units)
    scaling_factors = _prepare_scaling_factors(scaling_factors)
    data = _excavate_data(scenarios, variables, scaling_factors).to_dataframe().reset_index()
    data[["Scenario", "Variable", "Unit", "Value"]].to_csv(
        path_to_output,
        header=True,
        index=False
    )


def _excavate_data(scenarios, variables, scaling_factors):
    index = pd.MultiIndex.from_product(
        [
            scenarios,
            [variable.name for variable in variables]
        ],
        names=["Scenario", "Variable"]
    )
    data = pd.DataFrame(index=index, columns=["Unit", "Value"]).to_xarray()
    for variable in variables:
        if variable.value_function:
            _excavate_variable(data, scenarios, variable, scaling_factors)
    return data


def _excavate_variable(data, scenarios, variable, scaling_factors):
    values = {
        scenario_name: variable.value_function(scenario_data, scaling_factors)
        for scenario_name, scenario_data in scenarios.items()
    }
    for scenario, scenario_value in values.items():
        assert np.isscalar(scenario_value), scenario_value
        data.Value.loc[dict(Scenario=scenario, Variable=variable.name)] = scenario_value
        data.Unit.loc[dict(Scenario=scenario, Variable=variable.name)] = variable.unit


def _prepare_scaling_factors(scaling_factors):
    scaling_factors["energy"] = scaling_factors["power"] / 1e-6 # from MWh to TWh
    scaling_factors["stored_energy"] = scaling_factors["power"] / 1e-3 # from MWh to GWh
    scaling_factors["specific_cost"] = scaling_factors["monetary"] / scaling_factors["power"]
    scaling_factors["power"] = scaling_factors["power"] / 1e-3 # from MW to GW
    scaling_factors["monetary"] = scaling_factors["monetary"] / 1e-9 # from EUR to billion EUR
    return scaling_factors


def _set_up_variables(units):
    return [
        Variable("Cost|Total system", "billion EUR",
                 _cost_total_system),
        Variable("Cost|Capacity", "billion EUR",
                 _cost_capacity),
        Variable("Cost|Average total system", "EUR/MWh",
                 _levelised_total_cost),
        Variable("Cost|Variable", "billion EUR",
                 _cost_variable),
        Variable("Capacity|Wind onshore", "GW",
                 lambda model, sf: _capacity_for_tech(model, ONSHORE_WIND_TECHS, sf)),
        Variable("Capacity|Wind offshore", "GW",
                 lambda model, sf: _capacity_for_tech(model, ["wind_offshore"], sf)),
        Variable("Capacity|Wind", "GW",
                 lambda model, sf: _capacity_for_tech(model, ONSHORE_WIND_TECHS + ["wind_offshore"], sf)),
        Variable("Capacity|Solar PV", "GW",
                 lambda model, sf: _capacity_for_tech(model, PV_TECHS, sf)),
        Variable("Capacity|Hydro ROR", "GW",
                 lambda model, sf: _capacity_for_tech(model, ["hydro_run_of_river"], sf)),
        Variable("Capacity|Hydro Reservoir", "GW",
                 lambda model, sf: _capacity_for_tech(model, ["hydro_reservoir"], sf)),
        Variable("Capacity|Bioenergy", "GW",
                 lambda model, sf: _capacity_for_tech(model, ["biofuel"], sf)),
        Variable("Capacity|Storage|Short term|Power", "GW",
                 lambda model, sf: _capacity_for_tech(model, ["battery"], sf)),
        Variable("Capacity|Storage|Short term|Energy", "GWh",
                 lambda model, sf: _storage_capacity_for_tech(model, ["battery"], sf)),
        Variable("Capacity|Storage|Long term|Power", "GW",
                 lambda model, sf: _capacity_for_tech(model, ["hydrogen"], sf)),
        Variable("Capacity|Storage|Long term|Energy", "GWh",
                 lambda model, sf: _storage_capacity_for_tech(model, ["hydrogen"], sf)),
        Variable("Capacity|Storage|Pumped hydro|Power", "GW",
                 lambda model, sf: _capacity_for_tech(model, ["pumped_hydro"], sf)),
        Variable("Capacity|Storage|Pumped hydro|Energy", "GWh",
                 lambda model, sf: _storage_capacity_for_tech(model, ["pumped_hydro"], sf)),
        Variable("Capacity|Storage|Short and long term|Power", "GW",
                 lambda model, sf: _capacity_for_tech(model, ["battery", "hydrogen"], sf)),
        Variable("Capacity|Storage|Short and long term|Energy", "GWh",
                 lambda model, sf: _storage_capacity_for_tech(model, ["battery", "hydrogen"], sf)),
        Variable("Capacity|Gross import", "GW",
                 _transmission_capacity),
        Variable("Capacity|Gross export", "GW",
                 _transmission_capacity),
        Variable("Capacity|Transmission", "TW km",
                 _transmission_capacity_length),
        Variable("Energy|Wind onshore", "TWh",
                 lambda model, sf: _generation_for_tech(model, ONSHORE_WIND_TECHS, sf)),
        Variable("Energy|Wind offshore", "TWh",
                 lambda model, sf: _generation_for_tech(model, ["wind_offshore"], sf)),
        Variable("Energy|Solar PV", "TWh",
                 lambda model, sf: _generation_for_tech(model, PV_TECHS, sf)),
        Variable("Energy|Hydro ROR", "TWh",
                 lambda model, sf: _generation_for_tech(model, ["hydro_run_of_river"], sf)),
        Variable("Energy|Hydro Reservoir", "TWh",
                 lambda model, sf: _generation_for_tech(model, ["hydro_reservoir"], sf)),
        Variable("Energy|Bioenergy", "TWh",
                 lambda model, sf: _generation_for_tech(model, ["biofuel"], sf)),
        Variable("Energy|Renewable curtailment|Absolute", "TWh",
                 _absolute_curtailment),
        Variable("Energy|Renewable curtailment|Relative", "[%]",
                 _relative_curtailment),
        Variable("Energy|Storage|Short term", "TWh",
                 lambda model, sf: _generation_for_tech(model, ["battery"], sf)),
        Variable("Energy|Storage|Long term", "TWh",
                 lambda model, sf: _generation_for_tech(model, ["hydrogen"], sf)),
        Variable("Energy|Storage|Pumped hydro", "TWh",
                 lambda model, sf: _generation_for_tech(model, ["pumped_hydro"], sf)),
        Variable("Energy|Gross import", "TWh",
                 _transmission_generation),
        Variable("Energy|Gross export", "TWh",
                 _transmission_consumption),
        Variable("Energy|Gross import national level", "TWh",
                 lambda model, sf: _gross_import_national_level(model, units, sf)),
        Variable("Energy|Net import national level", "TWh",
                 lambda model, sf: _net_import_national_level(model, units, sf)),
    ]


@functools.lru_cache(maxsize=20, typed=False)
def _capacity(model, scaling_factor):
    return model.get_formatted_array("energy_cap") / scaling_factor


@functools.lru_cache(maxsize=CACHE_SIZE, typed=False)
def _storage_capacity(model, scaling_factor):
    return model.get_formatted_array("storage_cap") / scaling_factor


@functools.lru_cache(maxsize=CACHE_SIZE, typed=False)
def _generation(model, scaling_factor):
    return model.get_formatted_array("carrier_prod").squeeze("carriers") / scaling_factor


@functools.lru_cache(maxsize=CACHE_SIZE, typed=False)
def _consumption(model, scaling_factor):
    return model.get_formatted_array("carrier_con").squeeze("carriers") * (-1) / scaling_factor


@functools.lru_cache(maxsize=CACHE_SIZE, typed=False)
def _renewable_generation_potential(model, scaling_factor):
    resource = model.get_formatted_array("resource").sel(techs=RE_TECHS)
    capacity = model.get_formatted_array("energy_cap").sel(techs=RE_TECHS) / scaling_factor
    return (resource * capacity).sum(dim=["timesteps", "techs"])


@functools.lru_cache(maxsize=CACHE_SIZE, typed=False)
def _renewable_generation(model, scaling_factor):
    return (model.get_formatted_array("carrier_prod")
                 .sel(techs=RE_TECHS)
                 .sum(dim=["timesteps", "techs"])) / scaling_factor


def _cost_total_system(model, scaling_factors):
    return (model.get_formatted_array("cost")
                 .squeeze("costs")
                 .sum(["techs", "locs"])
                 .item()) / scaling_factors["monetary"]


def _cost_capacity(model, scaling_factors):
    return (model.get_formatted_array("cost_investment")
                 .squeeze("costs")
                 .sum(["techs", "locs"])
                 .item()) / scaling_factors["monetary"]


def _levelised_total_cost(model, scaling_factors):
    lcoe = (model.get_formatted_array("total_levelised_cost")
                 .squeeze(["carriers", "costs"])) / scaling_factors["specific_cost"]
    return lcoe.item()


def _cost_variable(model, scaling_factors):
    return (model.get_formatted_array("cost_var")
                 .squeeze("costs")
                 .sum(["techs", "timesteps", "locs"])
                 .item()) / scaling_factors["monetary"]


def _capacity_for_tech(model, techs, scaling_factors, optional=False):
    cap = _capacity(model, scaling_factors["power"])
    if optional and techs not in cap.techs:
        return 0
    else:
        return cap.sel(techs=techs).sum(["techs", "locs"]).item()


def _transmission_capacity(model, scaling_factors):
    capacity = _capacity(model, scaling_factors["power"])
    capacity = (capacity.sel(techs=capacity.techs.str.contains("ac_transmission"))
                        .sum(["techs", "locs"]))
    return capacity.item()


def _transmission_capacity_length(model, scaling_factors):
    capacity = _capacity(model, scaling_factors["power"]) * 1e-3 # to TW
    capacity = capacity.sel(techs=capacity.techs.str.contains("ac_transmission"))
    if "distance" in model._model_data:
        distance_km = model.get_formatted_array("distance") * M_TO_KM
        return (capacity * distance_km).sum(["techs", "locs"]).item() / 2 # counting each line twice
    else:
        assert capacity.sum() == 0
        return 0


def _storage_capacity_for_tech(model, techs, scaling_factors):
    return (_storage_capacity(model, scaling_factors["stored_energy"]).sel(techs=techs)
                                                                      .sum(["techs", "locs"])
                                                                      .item())


def _generation_for_tech(model, techs, scaling_factors, optional=False):
    generation = _generation(model, scaling_factors["energy"]).sum("timesteps")
    if optional and techs not in generation.techs:
        return 0
    return generation.sel(techs=techs).sum(["techs", "locs"]).item()


def _transmission_generation(model, scaling_factors):
    generation = _generation(model, scaling_factors["energy"]).sum("timesteps")
    return (generation.sel(techs=generation.techs.str.contains("ac_transmission"))
                      .sum(["techs", "locs"])
                      .item())


def _gross_import_national_level(model, units, scaling_factors):
    import_graph = _create_graph(model=model, units=units, scaling_factor=scaling_factors["energy"])

    def only_imports(node, neighbour):
        exchange = import_graph.adj[node][neighbour]["timeseries"][node][neighbour]
        return exchange.where(exchange > 0)

    return abs(sum([
        sum([only_imports(node, neighbour).sum("timesteps").item()
            for neighbour in import_graph.adj[node]])
        for node in import_graph.nodes
    ]))


def _net_import_national_level(model, units, scaling_factors):
    import_graph = _create_graph(model=model, units=units, scaling_factor=scaling_factors["energy"])
    national_exchanges = [
        sum([import_graph.adj[node][neighbour]["timeseries"][node][neighbour].sum("timesteps").item()
             for neighbour in import_graph.adj[node]])
        for node in import_graph.nodes
    ]
    return abs(sum([national_exchange for national_exchange in national_exchanges
                    if national_exchange >= 0]))


def _relative_curtailment(model, scaling_factors):
    potential = _renewable_generation_potential(model, scaling_factors["energy"]).sum("locs").item()
    generated = _renewable_generation(model, scaling_factors["energy"]).sum("locs").item()
    return (potential - generated) / potential * 100


def _absolute_curtailment(model, scaling_factors):
    potential = _renewable_generation_potential(model, scaling_factors["energy"]).sum("locs").item()
    generated = _renewable_generation(model, scaling_factors["energy"]).sum("locs").item()
    return potential - generated


def _transmission_consumption(model, scaling_factors):
    consumption = _consumption(model, scaling_factors["energy"]).sum("timesteps")
    return (consumption.sel(techs=consumption.techs.str.contains("ac_transmission"))
                       .sum(["techs", "locs"])
                       .item())


def _create_graph(model, units, scaling_factor): # FIXME almost duplicate from flows.py
    graph = nx.Graph()
    nodes = [country_code for country_code in units.country_code.unique()]
    graph.add_nodes_from(nodes)
    if "loc_techs_transmission" in model.inputs:
        edges = [(line.split(":")[0], line.split(":")[-1]) for line in model.inputs.loc_techs_transmission.values]
    else:
        edges = []
    edges = filter( # transform to linkes between countries
        lambda x: x[0] != x[1],
        map(
            lambda x: (units.loc[x[0], "country_code"], units.loc[x[1], "country_code"]),
            edges
        )
    )
    edges = list(set(tuple(sorted((a, b))) for a, b in edges)) # filter duplicates
    edges = [
        (source, sink, {
            "timeseries":
                {source:
                    {sink: _link_timeseries(model, units, source, sink, scaling_factor)},
                 sink:
                    {source: _link_timeseries(model, units, sink, source, scaling_factor)}
                 }
        })
        for (source, sink) in edges
    ]
    graph.add_edges_from(edges)
    return graph


def _link_timeseries(result, units, country_code_a, country_code_b, scaling_factor): # FIXME duplicate from flows.py
    gen = _transmission_carrier_prod(result, scaling_factor)
    con = _transmission_carrier_con(result, scaling_factor)
    nodes_a = units[units.country_code == country_code_a].index
    nodes_b = units[units.country_code == country_code_b].index
    a_to_b = con.sel(source=nodes_a, sink=nodes_b).sum(["sink", "source"])
    b_to_a = gen.sel(source=nodes_a, sink=nodes_b).sum(["sink", "source"])
    return a_to_b + b_to_a


@functools.lru_cache(maxsize=1, typed=False)
def _transmission_carrier_prod(model, scaling_factor): # FIXME almost duplicate from flows.py
    gen = _generation(model, scaling_factor)
    gen = gen.sel(techs=gen.techs.str.contains("ac_transmission"))
    return (gen.assign_coords(techs=[tech.item().split(":")[-1] for tech in gen.techs])
               .rename(locs="sink", techs="source"))


@functools.lru_cache(maxsize=1, typed=False)
def _transmission_carrier_con(model, scaling_factor): # FIXME almost duplicate from flows.py
    con = _consumption(model, scaling_factor) * -1
    con = con.sel(techs=con.techs.str.contains("ac_transmission"))
    return (con.assign_coords(techs=[tech.item().split(":")[-1] for tech in con.techs])
               .rename(locs="sink", techs="source"))


if __name__ == "__main__":
    excavate_all_results(
        paths_to_scenarios=snakemake.input.scenarios,
        path_to_units=snakemake.input.units,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
