import math

import calliope
import pandas as pd
import xarray as xr

EUR_PER_KW = 1 / 1e3 # from €/MW(h) to €/kW(h)
CT_PER_KW = 1e2 / 1e3 # from €/MW(h) to €ct/kW(h)
M_TO_1000KM = 1e-6
EPSILON = 1e-12
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
COST_SOURCES = {
    "open_field_pv": "[@JRC:2014] Table 7",
    "wind_onshore_monopoly": "[@JRC:2014] Table 4",
    "wind_offshore": "[@JRC:2014] Table 5",
    "biofuel": "[@JRC:2014] Table 48, [@RuizCastello:2015]",
    "ac_transmission": "[@JRC:2014] Table 39",
    "battery": "[@Schmidt:2019]",
    "hydrogen": "[@Schmidt:2019]",
    "load_shedding": ""
}


def main(path_to_model, scaling_factors, path_to_output):
    """Create table of important cost assumptions."""
    model = calliope.read_netcdf(path_to_model)
    eur_per_kw = scaling_factors["power"] / scaling_factors["monetary"] * EUR_PER_KW
    ct_per_kw = scaling_factors["power"] / scaling_factors["monetary"] * CT_PER_KW

    energy_cap = (model.get_formatted_array("cost_energy_cap")
                       .reindex(techs=list(TECHS.keys()))
                       .fillna(0)
                       .mean("locs")
                       .squeeze("costs")
                       .drop("costs")) * eur_per_kw
    energy_cap.loc["ac_transmission"] = transmission_cost(model, eur_per_kw)
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
        "Installation cost [€/kW]": energy_cap,
        "Installation cost [€/kWh]": storage_cap,
        "Annual cost [€/kW/yr]": annual_cost,
        "Variable cost [€ct/kWh]": variable_costs,
        "Lifetime [yr]": lifetime,
        "Source": pd.Series(COST_SOURCES).to_xarray().rename(index="techs")
    })
    all_costs.rename(techs="Technology").to_dataframe().rename(index=TECHS).to_csv(
        path_to_output,
        index=True,
        header=True,
        float_format="%.0f"
    )


def transmission_cost(model, scaling_factor):
    cost = model.get_formatted_array("cost_energy_cap").squeeze("costs") * scaling_factor
    distance = model.get_formatted_array("distance") * M_TO_1000KM
    rel_costs = (cost / distance).to_series().dropna()
    assert math.isclose(rel_costs.std(), 0, abs_tol=EPSILON)
    return rel_costs.iloc[0]


if __name__ == "__main__":
    main(
        path_to_model=snakemake.input.model,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
