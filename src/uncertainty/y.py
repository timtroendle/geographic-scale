import functools

import pandas as pd
import calliope

ENERGY_SCALING_FACTOR = 1e-3 # from MW(h) to GW(h)
PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
WIND_TECHS = ["wind_offshore", "wind_onshore_monopoly", "wind_onshore_competing"]
BIOFUEL_TECH = ["biofuel"]
STORAGE_TECHS = ["hydrogen", "battery"]


def determine_y(path_to_results, scaling_factors, experiment_id, path_to_output):
    results = calliope.read_netcdf(path_to_results)
    pd.DataFrame(
        data={
            "y-cost-eur": [system_cost(results, scaling_factors)],
            "y-pv-gw": [pv_capacity(results, scaling_factors)],
            "y-wind-gw": [wind_capacity(results, scaling_factors)],
            "y-biofuel-gw": [biofuel_capacity(results, scaling_factors)],
            "y-storage-gw": [storage_capacity_power(results, scaling_factors)],
            "y-storage-gwh": [storage_capacity_energy(results, scaling_factors)],
            "y-transmission-gwkm": [transmission_capacity(results, scaling_factors)],
        },
        index=[experiment_id]
    ).to_csv(path_to_output, sep="\t", index=True, header=True)


def system_cost(result, scaling_factors):
    return (result.results.cost
                  .squeeze("costs")
                  .sum("loc_techs_cost")
                  .item() / scaling_factors["monetary"])


def pv_capacity(result, scaling_factors):
    return (energy_cap(result, ENERGY_SCALING_FACTOR / scaling_factors["power"])
            .sel(techs=PV_TECHS)
            .sum(["techs", "locs"])
            .item())


def wind_capacity(result, scaling_factors):
    return (energy_cap(result, ENERGY_SCALING_FACTOR / scaling_factors["power"])
            .sel(techs=WIND_TECHS)
            .sum(["techs", "locs"])
            .item())


def biofuel_capacity(result, scaling_factors):
    return (energy_cap(result, ENERGY_SCALING_FACTOR / scaling_factors["power"])
            .sel(techs=BIOFUEL_TECH)
            .sum(["locs"])
            .item())


def storage_capacity_power(result, scaling_factors):
    return (energy_cap(result, ENERGY_SCALING_FACTOR / scaling_factors["power"])
            .sel(techs=STORAGE_TECHS)
            .sum(["techs", "locs"])
            .item())


def storage_capacity_energy(result, scaling_factors):
    return (result.get_formatted_array("storage_cap")
                  .sel(techs=STORAGE_TECHS)
                  .sum(["techs", "locs"])
                  .item() / scaling_factors["power"] * ENERGY_SCALING_FACTOR)


def transmission_capacity(result, scaling_factors):
    energy_cap_gw = energy_cap(result, ENERGY_SCALING_FACTOR / scaling_factors["power"])
    trans_cap_gw = energy_cap_gw.sel(techs=["ac_transmission" in str(tech) for tech in energy_cap_gw.techs])
    if "distance" in result._model_data:
        distance_km = result.get_formatted_array("distance") / 1000
        return (trans_cap_gw * distance_km).sum(["techs", "locs"]).item() / 2 # counting each line twice
    else:
        assert trans_cap_gw.sum() == 0
        return 0


@functools.lru_cache(maxsize=1, typed=False)
def energy_cap(model, scaling_factor):
    return model.get_formatted_array("energy_cap") * scaling_factor


if __name__ == "__main__":
    determine_y(
        path_to_results=snakemake.input.result,
        scaling_factors=snakemake.params.scaling_factors,
        experiment_id=snakemake.params.experiment_id,
        path_to_output=snakemake.output[0]
    )
