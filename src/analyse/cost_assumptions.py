import calliope
import xarray as xr

OUTPUT_SCALING_FACTOR = 1 / 1e3 # from €/MW(h) to €/kW(h)
TECHS = {
    "open_field_pv": "PV",
    "wind_onshore_monopoly": "Onshore wind",
    "wind_offshore": "Offshore wind",
    "ac_transmission": "AC transmission^",
    "battery": "Battery",
    "hydrogen": "Hydrogen",
    "load_shedding": "Load shedding"
}


def main(path_to_model, scaling_factors, path_to_output):
    """Create table of important cost assumptions."""
    # TODO add ac transmission installation costs
    model = calliope.read_netcdf(path_to_model)
    cost_scaling_factor = scaling_factors["power"] / scaling_factors["monetary"] * OUTPUT_SCALING_FACTOR

    energy_cap = (model.get_formatted_array("cost_energy_cap")
                       .reindex(techs=list(TECHS.keys()))
                       .fillna(0)
                       .mean("locs")
                       .squeeze("costs")
                       .drop("costs")) * cost_scaling_factor
    annual_cost = (model.get_formatted_array("cost_om_annual")
                        .reindex(techs=list(TECHS.keys()))
                        .fillna(0)
                        .mean("locs")
                        .squeeze("costs")
                        .drop("costs")) * cost_scaling_factor
    storage_cap = (model.get_formatted_array("cost_storage_cap")
                        .reindex(techs=list(TECHS.keys()))
                        .fillna(0)
                        .mean("locs")
                        .squeeze("costs")
                        .drop("costs")) * cost_scaling_factor
    lifetime = (model.get_formatted_array("lifetime")
                     .reindex(techs=list(TECHS.keys()))
                     .fillna(0)
                     .mean("locs"))
    variable_costs = (model.get_formatted_array("cost_om_prod")
                           .reindex(techs=list(TECHS.keys()))
                           .fillna(0)
                           .mean("locs")
                           .squeeze("costs")
                           .drop("costs")) * cost_scaling_factor

    all_costs = xr.Dataset({
        "installation cost [€/kW]": energy_cap,
        "installation cost [€/kWh]": storage_cap,
        "annual cost [EUR/kW/yr]": annual_cost,
        "variable cost [EUR/kWh]": variable_costs,
        "lifetime [yr]": lifetime
    })
    all_costs.rename(techs="technology").to_dataframe().rename(index=TECHS).to_csv(
        path_to_output,
        index=True,
        header=True,
        float_format="%.0f"
    )


if __name__ == "__main__":
    main(
        path_to_model=snakemake.input.model,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
