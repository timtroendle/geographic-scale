from collections import OrderedDict

import pandas as pd

AUTARKY_SCALE = '<i class="fas fa-layer-group"></i>'
AUTARKY_LEVEL = '<i class="fas fa-shield-alt"></i>'
GRID_SCALE = '<i class="fab fa-connectdevelop"></i>'
SOLAR = '<i class="fas fa-sun"></i>'
WIND = '<i class="fas fa-wind"></i>'
BIOFUEL = '<i class="fas fa-leaf"></i>'
STORAGE = '<i class="fas fa-battery-three-quarters"></i>'
TRANSMISSION = '<i class="fab fa-connectdevelop"></i>'
IMPORT = '<i class="fas fa-shopping-cart"></i>'
CURTAILMENT = '<i class="fas fa-traffic-light"></i>'
LOAD_SHED = '<i class="fas fa-truck-loading"></i>'

VARIABLE_MAP = OrderedDict([
    ("Capacity|Solar PV", f"{SOLAR} [GW]"),
    ("Capacity|Wind", f"{WIND} [GW]"),
    ("Capacity|Bioenergy", f"{BIOFUEL} [GW]"),
    ("Capacity|Storage|Short and long term|Power", f"{STORAGE} [GW]"),
    ("Capacity|Storage|Short and long term|Energy", f"{STORAGE} [GWh]"),
    ("Capacity|Transmission", f"{TRANSMISSION} [TW km]"),
    ("Energy|Gross import national level", f"{IMPORT} gross [TWh]"),
    ("Energy|Net import national level", f"{IMPORT} net [TWh]"),
    ("Energy|Renewable curtailment|Relative", f"{CURTAILMENT} [%]"),
    ("Energy|Load shedding|Relative", f"{LOAD_SHED} [‰]")
])


def main(path_to_aggregated_results, path_to_output):
    results = pd.read_csv(path_to_aggregated_results, index_col=[0, 1]).to_xarray()
    results = results.sel(Variable=list(VARIABLE_MAP.keys()))
    results = (results.assign_coords(Variable=[VARIABLE_MAP[old] for old in results.Variable.values])
                      .Value
                      .to_dataframe()
                      .reset_index()
                      .pivot(index="Scenario", columns="Variable", values="Value")[list(VARIABLE_MAP.values())]
                      .rename(index=nice_scenario_name)
                      .to_csv(path_to_output, index=True, header=True, float_format="%.0f"))


def nice_scenario_name(scenario_name):
    autarky_scale, _, autarky_level, grid_scale, _ = scenario_name.split("-")
    return f"""{AUTARKY_SCALE}: {autarky_scale}<br>
    {AUTARKY_LEVEL}: {autarky_level}%<br>
    {GRID_SCALE}: {grid_scale}<br>
    """


if __name__ == "__main__":
    main(
        path_to_aggregated_results=snakemake.input.results,
        path_to_output=snakemake.output[0]
    )
