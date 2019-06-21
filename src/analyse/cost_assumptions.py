import calliope
import xarray as xr

EUR_PER_KW = 1 / 1e3 # from €/MW(h) to €/kW(h)
CT_PER_KW = 1e2 / 1e3 # from €/MW(h) to €ct/kW(h)
TECHS = {
    "open_field_pv": "PV",
    "wind_onshore_monopoly": "Onshore wind",
    "wind_offshore": "Offshore wind",
    "biofuel": "Biofuel",
    "ac_transmission": "AC transmission^",
    "battery": "Short term storage",
    "hydrogen": "Long term storage",
    "load_shedding": "Load shedding"
}


def main(path_to_model, scaling_factors, path_to_output):
    """Create table of important cost assumptions."""
    # TODO add ac transmission installation costs
    model = calliope.read_netcdf(path_to_model)
    eur_per_kw = scaling_factors["power"] / scaling_factors["monetary"] * EUR_PER_KW
    ct_per_kw = scaling_factors["power"] / scaling_factors["monetary"] * CT_PER_KW

    energy_cap = (model.get_formatted_array("cost_energy_cap")
                       .reindex(techs=list(TECHS.keys()))
                       .fillna(0)
                       .mean("locs")
                       .squeeze("costs")
                       .drop("costs")) * eur_per_kw
    annual_cost = (model.get_formatted_array("cost_om_annual")
                        .reindex(techs=list(TECHS.keys()))
                        .fillna(0)
                        .mean("locs")
                        .squeeze("costs")
                        .drop("costs")) * eur_per_kw
    storage_cap = (model.get_formatted_array("cost_storage_cap")
                        .reindex(techs=list(TECHS.keys()))
                        .fillna(0)
                        .mean("locs")
                        .squeeze("costs")
                        .drop("costs")) * eur_per_kw
    lifetime = (model.get_formatted_array("lifetime")
                     .reindex(techs=list(TECHS.keys()))
                     .fillna(0)
                     .mean("locs"))
    variable_costs_prod = (model.get_formatted_array("cost_om_prod")
                                .reindex(techs=list(TECHS.keys()))
                                .fillna(0)
                                .mean("locs")
                                .squeeze("costs")
                                .drop("costs")) * ct_per_kw
    variable_costs_con = (model.get_formatted_array("cost_om_con")
                               .reindex(techs=list(TECHS.keys()))
                               .fillna(0)
                               .mean("locs")
                               .squeeze("costs")
                               .drop("costs")) * ct_per_kw
    variable_costs = variable_costs_prod + variable_costs_con

    all_costs = xr.Dataset({
        "installation cost [€/kW]": energy_cap,
        "installation cost [€/kWh]": storage_cap,
        "annual cost [€/kW/yr]": annual_cost,
        "variable cost [€ct/kWh]": variable_costs,
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
