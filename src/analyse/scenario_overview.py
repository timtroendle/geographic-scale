from collections import OrderedDict

import pandas as pd

SOLAR = 'PV'
WIND = 'Wind'
BIOFUEL = 'Bioenergy'
SHORT_TERM_STORAGE = 'Battery'
LONG_TERM_STORAGE = 'Hydrogen'
TRANSMISSION = 'Transmission'
IMPORT = 'National import'
CURTAILMENT = 'Curtailment'

SCENARIOS = [
    "continental-autarky-100-continental-grid",
    "national-autarky-100-national-grid",
    "regional-autarky-100-regional-grid",
    "national-autarky-100-continental-grid",
    "national-autarky-85-continental-grid",
    "national-autarky-70-continental-grid",
    "regional-autarky-100-continental-grid",
    "regional-autarky-85-continental-grid",
    "regional-autarky-70-continental-grid",
    "regional-autarky-100-national-grid",
    "regional-autarky-85-national-grid",
    "regional-autarky-70-national-grid",

]
VARIABLES_TABLE_1 = OrderedDict([
    ("Capacity|Solar PV", f"{SOLAR} (GW)"),
    ("Capacity|Wind", f"{WIND} (GW)"),
    ("Energy|Renewable curtailment|Relative total", f"{CURTAILMENT} (%)"),
    ("Capacity|Bioenergy", f"{BIOFUEL} (GW)"),
    ("Capacity|Storage|Short term|Power", f"{SHORT_TERM_STORAGE} (GW)"),
    ("Capacity|Storage|Short term|Energy", f"{SHORT_TERM_STORAGE} (GWh)"),
    ("Capacity|Storage|Long term|Power", f"{LONG_TERM_STORAGE} (GW)"),
    ("Capacity|Storage|Long term|Energy", f"{LONG_TERM_STORAGE} (GWh)"),
])
VARIABLES_TABLE_2 = OrderedDict([
    ("Capacity|Transmission", f"{TRANSMISSION} (TW km)"),
    ("Energy|Gross import national level", f"{IMPORT} gross (TWh)"),
    ("Energy|Net import national level", f"{IMPORT} net (TWh)")
])


def main(path_to_aggregated_results, path_to_output_table1, path_to_output_table2):
    results = pd.read_csv(path_to_aggregated_results, index_col=[0, 1]).to_xarray()
    all_variables = [VARIABLES_TABLE_1, VARIABLES_TABLE_2]
    all_output_paths = [path_to_output_table1, path_to_output_table2]
    for variables, path_to_output in zip(all_variables, all_output_paths):
        filtered_results = results.sel(Variable=list(variables.keys()))
        (
            filtered_results
            .assign_coords(Variable=[variables[old] for old in filtered_results.Variable.values])
            .sel(Scenario=SCENARIOS)
            .Value
            .to_dataframe()
            .reset_index()
            .pivot(index="Scenario", columns="Variable", values="Value")
            .loc[SCENARIOS, list(variables.values())]
            .rename(index=nice_scenario_name)
            .rename_axis(index="Case")
            .to_csv(path_to_output, index=True, header=True, float_format="%.0f")
        )


def nice_scenario_name(scenario_name):
    autarky_scale, _, autarky_level, grid_scale, _ = scenario_name.split("-")
    if autarky_scale == grid_scale:
        return f"{autarky_scale.capitalize()} base case"
    else:
        return (f"{grid_scale.capitalize()} scale with {autarky_scale} self-sufficiency "
                f"and {100 - int(autarky_level)}% net imports")


if __name__ == "__main__":
    main(
        path_to_aggregated_results=snakemake.input.results,
        path_to_output_table1=snakemake.output.table1,
        path_to_output_table2=snakemake.output.table2
    )
