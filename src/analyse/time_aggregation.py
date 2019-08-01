"""Collect time aggregated results."""
from pathlib import Path

import calliope
import geopandas as gpd
import xarray as xr

VARIABLE_SCALING_FACTOR = {
    "energy_cap": lambda sf: 1 / sf["power"],
    "storage_cap": lambda sf: 1 / sf["power"],
    "carrier_prod": lambda sf: 1 / sf["power"],
    "carrier_con": lambda sf: 1 / sf["power"],
    "total_levelised_cost": lambda sf: sf["power"] / sf["monetary"],
    "capacity_factor": lambda sf: 1,
    "resource": lambda sf: 1,
    "cost": lambda sf: 1 / sf["monetary"]
}


def excavate_all_results(paths_to_scenarios, path_to_units, scaling_factors, path_to_output):
    """Collect time aggregated results of all scenarios."""
    scenarios = {
        Path(path_to_scenario).parent.name: calliope.read_netcdf(path_to_scenario)
        for path_to_scenario in paths_to_scenarios
    }
    units = (gpd.read_file(path_to_units)
                .set_index("id")
                .rename(index=lambda idx: idx.replace(".", "-"))
                .rename_axis(index="locs")
                .to_xarray())
    exporter = CalliopeExporter(scenarios, units, scaling_factors)
    ds = xr.Dataset({
        variable_name: exporter(variable_name)
        for variable_name in VARIABLE_SCALING_FACTOR.keys()
    })
    ds.coords["country_code"] = units.country_code
    ds.coords["tech_group"] = list(scenarios.values())[0].get_formatted_array("inheritance")
    ds.coords["tech_group"].loc[:] = [parent.split('.')[-1] for parent in ds["tech_group"].values]
    ds.to_netcdf(path_to_output)


class CalliopeExporter:

    def __init__(self, scenarios, units, scaling_factors):
        self.__scenarios = scenarios
        self.__units = units
        self.__scaling_factors = scaling_factors

    def __call__(self, variable_name):
        if variable_name not in VARIABLE_SCALING_FACTOR.keys():
            raise ValueError(f"Unknown variable name: {variable_name}.")
        return xr.concat([
            self._read(model, variable_name).expand_dims(scenario=[name], axis=0)
            for name, model in self.__scenarios.items()
        ], dim="scenario")

    def _read(self, model, variable_name):
        data = model.get_formatted_array(variable_name)
        if "carriers" in data.dims:
            data = data.squeeze("carriers").drop("carriers")
        if "costs" in data.dims:
            data = data.squeeze("costs").drop("costs")
        if "timesteps" in data.dims:
            data = data.sum("timesteps")
        if "techs" in data.dims:
            if data.techs.str.contains("ac_transmission").any():
                data = (data
                        .groupby(data.techs.where(~data.techs.str.contains("ac_transmission"), "ac_transmission"))
                        .sum("techs"))
            data = data.dropna(dim="techs", how="all")
        return data * VARIABLE_SCALING_FACTOR[variable_name](self.__scaling_factors)


if __name__ == "__main__":
    excavate_all_results(
        paths_to_scenarios=snakemake.input.scenarios,
        path_to_units=snakemake.input.units,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
