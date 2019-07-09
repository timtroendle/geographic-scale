import functools
from pathlib import Path

import calliope
import pandas as pd
import geopandas as gpd
import networkx as nx

ENERGY_SCALING_FACTOR = 1 / 1e3 # from MW(h) to GW(h)
MWH_TO_TWH = 1e-6
CACHE_SIZE = 20
PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
WIND_TECHS = ["wind_offshore", "wind_onshore_monopoly", "wind_onshore_competing"]
HYDRO_TECHS = ["hydro_run_of_river"]
BIOMASS_TECHS = ["biofuel"]
CURTAILABLE_RE_TECHS = PV_TECHS + WIND_TECHS + HYDRO_TECHS
ALL_STORAGE_TECHS = ["battery", "biofuel", "hydro_reservoir", "hydrogen", "pumped_hydro"]
DEPLOYABLE_STORAGE_TECHS = ["battery", "hydrogen"]
DEMAND_TECH = "demand_elec"
AUTARKY_SCALE = '<i class="fas fa-layer-group"></i>'
AUTARKY_LEVEL = '<i class="fas fa-shield-alt"></i>'
GRID_SCALE = '<i class="fab fa-connectdevelop"></i>'


def main(paths_to_results, path_to_units, scaling_factors, path_to_output):
    scenario_results = {
        scenario_name(path_to_results): calliope.read_netcdf(path_to_results)
        for path_to_results in paths_to_results
    }
    units = (gpd.read_file(path_to_units)
                .set_index("id")
                .rename(index=lambda idx: idx.replace(".", "-"))
                .rename_axis(index="locs"))
    import_graphs = [
        create_graph(result, units, scaling_factors["power"])
        for result in scenario_results.values()
    ]
    data = pd.DataFrame(
        index=scenario_results.keys(),
        data={
            "PV capacity [GW]": [pv_capacity(result, scaling_factors)
                                 for result in scenario_results.values()],
            "Wind capacity [GW]": [wind_capacity(result, scaling_factors)
                                   for result in scenario_results.values()],
            "Biomass capacity [GW]": [biomass_capacity(result, scaling_factors)
                                      for result in scenario_results.values()],
            "Storage capacity [GW]": [storage_capacity_power(result, scaling_factors)
                                      for result in scenario_results.values()],
            "Storage capacity [GWh]": [storage_capacity_energy(result, scaling_factors)
                                       for result in scenario_results.values()],
            "Transmission capacity [GW km]": [transmission_capacity(result, scaling_factors)
                                              for result in scenario_results.values()],
            "Gross imports [TWh]": [gross_import(import_graph)
                                    for import_graph in import_graphs],
            "Net imports [TWh]": [net_import(import_graph)
                                  for import_graph in import_graphs],
            "Curtailment [%]": [relative_curtailment(result, scaling_factors)
                                for result in scenario_results.values()],
            "Load shed [‰]": [shed_load(result, scaling_factors)
                              for result in scenario_results.values()]
        }
    )
    data.index.name = "Scenario"
    data.to_csv(
        path_to_output,
        index=True,
        header=True,
        float_format="%.0f"
    )


def scenario_name(path_to_result):
    autarky_scale, _, autarky_level, grid_scale, _ = Path(path_to_result).parent.name.split("-")
    return f"""{AUTARKY_SCALE}: {autarky_scale}<br>
    {AUTARKY_LEVEL}: {autarky_level}%<br>
    {GRID_SCALE}: {grid_scale}<br>
    """


def pv_capacity(result, scaling_factors):
    return capacity(result, PV_TECHS, scaling_factors)


def wind_capacity(result, scaling_factors):
    return capacity(result, WIND_TECHS, scaling_factors)


def biomass_capacity(result, scaling_factors):
    return capacity(result, BIOMASS_TECHS, scaling_factors)


def storage_capacity_power(result, scaling_factors):
    return capacity(result, DEPLOYABLE_STORAGE_TECHS, scaling_factors)


def capacity(result, techs, scaling_factors):
    return (result.get_formatted_array("energy_cap")
                  .sel(techs=techs)
                  .sum()
                  .item() / scaling_factors["power"] * ENERGY_SCALING_FACTOR)


def storage_capacity_energy(result, scaling_factors):
    storage_cap = result.get_formatted_array("storage_cap")
    assert set(storage_cap.techs.values) == set(ALL_STORAGE_TECHS)
    return (storage_cap.sel(techs=DEPLOYABLE_STORAGE_TECHS)
                       .sum(["techs", "locs"])
                       .item() / scaling_factors["power"] * ENERGY_SCALING_FACTOR)


def transmission_capacity(result, scaling_factors):
    energy_cap = result.get_formatted_array("energy_cap") / scaling_factors["power"]
    trans_cap_mw = energy_cap.sel(techs=["ac_transmission" in str(tech) for tech in energy_cap.techs])
    if "distance" in result._model_data:
        distance_km = result.get_formatted_array("distance") / 1000
        return (trans_cap_mw * distance_km).sum().item() * ENERGY_SCALING_FACTOR / 2 # counting each line twice
    else:
        assert trans_cap_mw.sum() == 0
        return 0


@functools.lru_cache(maxsize=CACHE_SIZE, typed=False)
def _generation(result, scaling_factor):
    return result.get_formatted_array("carrier_prod").squeeze("carriers") / scaling_factor


def gross_import(import_graph):
    return sum([
        abs(import_graph.edges[edge]["timeseries"][edge[0]][edge[1]]).sum("timesteps").item()
        for edge in import_graph.edges
    ])


def net_import(import_graph):
    return sum([
        abs(import_graph.edges[edge]["timeseries"][edge[0]][edge[1]].sum("timesteps").item())
        for edge in import_graph.edges
    ])


def relative_curtailment(result, scaling_factors):
    resource = result.get_formatted_array("resource").sel(techs=CURTAILABLE_RE_TECHS)
    capacity = result.get_formatted_array("energy_cap").sel(techs=CURTAILABLE_RE_TECHS) / scaling_factors["power"]
    potential = (resource * capacity).sum(dim=["timesteps", "techs", "locs"])
    generated = (_generation(result, scaling_factors["power"]).sel(techs=CURTAILABLE_RE_TECHS)
                                                              .sum(dim=["timesteps", "techs", "locs"]))
    return ((potential - generated) / potential).item() * 100


def shed_load(result, scaling_factors):
    gen = _generation(result, scaling_factors["power"])
    if "load_shedding" in gen.techs:
        shed = gen.sel(techs="load_shedding").sum().item()
        demand = (result.get_formatted_array("carrier_con")
                        .sel(techs=DEMAND_TECH)
                        .sum()
                        .item()) * -1 * scaling_factors["power"]
        return (shed / demand * 1000)
    else:
        return 0


def create_graph(model, units, scaling_factor): # FIXME almost duplicate from flows.py
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
                    {sink: link_timeseries(model, units, source, sink, scaling_factor)},
                 sink:
                    {source: link_timeseries(model, units, sink, source, scaling_factor)}
                 }
        })
        for (source, sink) in edges
    ]
    graph.add_edges_from(edges)
    return graph


def link_timeseries(result, units, country_code_a, country_code_b, scaling_factor): # FIXME duplicate from flows.py
    gen = transmission_carrier_prod(result, scaling_factor)
    nodes_a = units[units.country_code == country_code_a].index
    nodes_b = units[units.country_code == country_code_b].index
    a_to_b = gen.sel(source=nodes_a, sink=nodes_b).sum(["sink", "source"])
    b_to_a = gen.sel(source=nodes_b, sink=nodes_a).sum(["sink", "source"])
    return a_to_b - b_to_a


@functools.lru_cache(maxsize=1, typed=False)
def transmission_carrier_prod(model, scaling_factor): # FIXME almost duplicate from flows.py
    gen = _generation(model, scaling_factor) * MWH_TO_TWH
    gen = gen.sel(techs=gen.techs.str.contains("ac_transmission"))
    return (gen.assign_coords(techs=[tech.item().split(":")[-1] for tech in gen.techs])
               .rename(locs="sink", techs="source"))


if __name__ == "__main__":
    main(
        paths_to_results=snakemake.input.results,
        path_to_units=snakemake.input.units,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
