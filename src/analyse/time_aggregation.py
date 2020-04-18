"""Collect time aggregated results."""
from dataclasses import dataclass

import calliope
import geopandas as gpd
import xarray as xr


@dataclass
class Variable:
    name: str
    scaling_factor: callable
    unit: str
    description: str


VARIABLES = [
    Variable(
        name="energy_cap",
        scaling_factor=lambda sf: 1 / sf["power"],
        unit="MW",
        description="installed capacity"
    ),
    Variable(
        name="storage_cap",
        scaling_factor=lambda sf: 1 / sf["power"],
        unit="MWh",
        description="installed storage capacity"
    ),
    Variable(
        name="carrier_prod",
        scaling_factor=lambda sf: 1 / sf["power"],
        unit="MWh",
        description="generated electricity"
    ),
    Variable(
        name="carrier_con",
        scaling_factor=lambda sf: 1 / sf["power"],
        unit="MWh",
        description="consumed electricity"
    ),
    Variable(
        name="total_levelised_cost",
        scaling_factor=lambda sf: sf["power"] / sf["monetary"],
        unit="EUR/MWh",
        description="total levelised cost"
    ),
    Variable(
        name="capacity_factor",
        scaling_factor=lambda sf: 1,
        unit="-",
        description="capacity factors of supply, storage, and transmission"
    ),
    Variable(
        name="resource",
        scaling_factor=lambda sf: 1,
        unit="-",
        description="resource that technology is based upon"
    ),
    Variable(
        name="cost",
        scaling_factor=lambda sf: 1 / sf["monetary"],
        unit="EUR",
        description="technology cost"
    )
]


def excavate_all_results(paths_to_scenarios, path_to_units, scaling_factors, path_to_output):
    """Collect time aggregated results of all scenarios."""
    scenarios = {
        calliope.read_netcdf(path_to_scenario).results.attrs["scenario"]: calliope.read_netcdf(path_to_scenario)
        for path_to_scenario in paths_to_scenarios
    }
    units = (gpd.read_file(path_to_units)
                .set_index("id")
                .rename(index=lambda idx: idx.replace(".", "-"))
                .rename_axis(index="locs")
                .to_xarray())
    exporter = CalliopeExporter(scenarios, units, scaling_factors)
    ds = xr.Dataset({
        variable.name: exporter(variable)
        for variable in VARIABLES
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

    def __call__(self, variable):
        if variable not in VARIABLES:
            raise ValueError(f"Unknown variable: {variable}.")
        return xr.concat([
            self._read(model, variable).expand_dims(scenario=[name], axis=0)
            for name, model in self.__scenarios.items()
        ], dim="scenario")

    def _read(self, model, variable):
        data = model.get_formatted_array(variable.name)
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
        data = data * variable.scaling_factor(self.__scaling_factors)
        data.attrs["unit"] = variable.unit
        data.attrs["description"] = variable.description
        return data


if __name__ == "__main__":
    excavate_all_results(
        paths_to_scenarios=snakemake.input.scenarios,
        path_to_units=snakemake.input.units,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
