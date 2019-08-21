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

SCENARIOS = [
    "continental-autarky-100-continental-grid",
    "national-autarky-100-continental-grid",
    "national-autarky-85-continental-grid",
    "national-autarky-70-continental-grid",
    "national-autarky-100-national-grid",
    "regional-autarky-100-continental-grid",
    "regional-autarky-85-continental-grid",
    "regional-autarky-70-continental-grid",
    "regional-autarky-100-national-grid",
    "regional-autarky-85-national-grid",
    "regional-autarky-70-national-grid",
    "regional-autarky-100-regional-grid",

]
VARIABLES_TABLE_1 = OrderedDict([
    ("Capacity|Solar PV", f"{SOLAR} [GW]"),
    ("Capacity|Wind", f"{WIND} [GW]"),
    ("Energy|Renewable curtailment|Relative total", f"{CURTAILMENT} [%]"),
    ("Capacity|Bioenergy", f"{BIOFUEL} [GW]"),
    ("Capacity|Storage|Short term|Power", f"{STORAGE} short [GW]"),
    ("Capacity|Storage|Short term|Energy", f"{STORAGE} short [GWh]"),
    ("Capacity|Storage|Long term|Power", f"{STORAGE} long [GW]"),
    ("Capacity|Storage|Long term|Energy", f"{STORAGE} long [GWh]"),
])
VARIABLES_TABLE_2 = OrderedDict([
    ("Capacity|Transmission", f"{TRANSMISSION} [TW km]"),
    ("Energy|Gross import national level", f"{IMPORT} gross [TWh]"),
    ("Energy|Net import national level", f"{IMPORT} net [TWh]")
])


def main(path_to_aggregated_results, path_to_output_table1, path_to_output_table2):
    results = pd.read_csv(path_to_aggregated_results, index_col=[0, 1]).to_xarray()
    results_1 = results.sel(Variable=list(VARIABLES_TABLE_1.keys()))
    (results_1.assign_coords(Variable=[VARIABLES_TABLE_1[old] for old in results_1.Variable.values])
              .sel(Scenario=SCENARIOS)
              .Value
              .to_dataframe()
              .reset_index()
              .pivot(index="Scenario", columns="Variable", values="Value")
              .loc[SCENARIOS, list(VARIABLES_TABLE_1.values())]
              .rename(index=nice_scenario_name)
              .rename_axis(index="Layout")
              .to_csv(path_to_output_table1, index=True, header=True, float_format="%.0f"))
    results_2 = results.sel(Variable=list(VARIABLES_TABLE_2.keys()))
    (results_2.assign_coords(Variable=[VARIABLES_TABLE_2[old] for old in results_2.Variable.values])
              .sel(Scenario=SCENARIOS)
              .Value
              .to_dataframe()
              .reset_index()
              .pivot(index="Scenario", columns="Variable", values="Value")
              .loc[SCENARIOS, list(VARIABLES_TABLE_2.values())]
              .rename(index=nice_scenario_name)
              .rename_axis(index="Layout")
              .to_csv(path_to_output_table2, index=True, header=True, float_format="%.0f"))


def nice_scenario_name(scenario_name):
    autarky_scale, _, autarky_level, grid_scale, _ = scenario_name.split("-")
    return f"""{AUTARKY_SCALE}: {autarky_scale.capitalize()}<br>
    {AUTARKY_LEVEL}: {100 - int(autarky_level)}%<br>
    {GRID_SCALE}: {grid_scale.capitalize()}<br>
    """


if __name__ == "__main__":
    main(
        path_to_aggregated_results=snakemake.input.results,
        path_to_output_table1=snakemake.output.table1,
        path_to_output_table2=snakemake.output.table2
    )
