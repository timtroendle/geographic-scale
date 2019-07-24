import functools

import pandas as pd
import calliope

ENERGY_SCALING_FACTOR = 1e-3 # from MW(h) to GW(h)
PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
WIND_TECHS = ["wind_offshore", "wind_onshore_monopoly", "wind_onshore_competing"]
BIOFUEL_TECH = ["biofuel"]
STORAGE_TECHS = ["hydrogen", "battery"]


def determine_y(path_to_results_of_large_scale, path_to_results_of_small_scale,
                scaling_factors, experiment_id, path_to_output):
    large_scale = calliope.read_netcdf(path_to_results_of_large_scale)
    small_scale = calliope.read_netcdf(path_to_results_of_small_scale)
    costs_large_scale = system_cost(large_scale, scaling_factors)
    costs_small_scale = system_cost(small_scale, scaling_factors)
    cost_diff = costs_small_scale - costs_large_scale
    cost_diff_relative = (costs_small_scale - costs_large_scale) / costs_large_scale
    pd.DataFrame(
        data={
            "y-large-scale-cost-eur": [costs_large_scale],
            "y-small-scale-cost-eur": [costs_small_scale],
            "y-cost-diff-eur": [cost_diff],
            "y-cost-diff-relative": [cost_diff_relative],
            "y-large-scale-pv-gw": [pv_capacity(large_scale, scaling_factors)],
            "y-small-scale-pv-gw": [pv_capacity(small_scale, scaling_factors)],
            "y-large-scale-wind-gw": [wind_capacity(large_scale, scaling_factors)],
            "y-small-scale-wind-gw": [wind_capacity(small_scale, scaling_factors)],
            "y-large-scale-biofuel-gw": [biofuel_capacity(large_scale, scaling_factors)],
            "y-small-scale-biofuel-gw": [biofuel_capacity(small_scale, scaling_factors)],
            "y-large-scale-storage-gw": [storage_capacity_power(large_scale, scaling_factors)],
            "y-small-scale-storage-gw": [storage_capacity_power(small_scale, scaling_factors)],
            "y-large-scale-storage-gwh": [storage_capacity_energy(large_scale, scaling_factors)],
            "y-small-scale-storage-gwh": [storage_capacity_energy(small_scale, scaling_factors)],
            "y-large-scale-transmission-gwkm": [transmission_capacity(large_scale, scaling_factors)],
            "y-small-scale-transmission-gwkm": [transmission_capacity(small_scale, scaling_factors)]
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
        path_to_results_of_large_scale=snakemake.input.large_scale,
        path_to_results_of_small_scale=snakemake.input.small_scale,
        scaling_factors=snakemake.params.scaling_factors,
        experiment_id=snakemake.params.experiment_id,
        path_to_output=snakemake.output[0]
    )
