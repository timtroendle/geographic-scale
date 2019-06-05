import pandas as pd
import calliope

ENERGY_SCALING_FACTOR = 1 / 1e3 # from MW(h) to GW(h)
PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
WIND_TECHS = ["wind_offshore", "wind_onshore_monopoly", "wind_onshore_competing"]


def determine_y(path_to_results_of_large_scale, path_to_results_of_small_scale,
                scaling_factors, experiment_id, path_to_output):
    large_scale = calliope.read_netcdf(path_to_results_of_large_scale)
    small_scale = calliope.read_netcdf(path_to_results_of_small_scale)
    costs_large_scale = system_cost(large_scale, scaling_factors)
    costs_small_scale = system_cost(small_scale, scaling_factors)
    cost_diff = costs_small_scale - costs_large_scale
    pd.DataFrame(
        data={
            "y-large-scale-cost-eur": [costs_large_scale],
            "y-small-scale-cost-eur": [costs_small_scale],
            "y-cost-diff-eur": [cost_diff],
            "y-large-scale-pv-gw": [pv_capacity(large_scale, scaling_factors)],
            "y-small-scale-pv-gw": [pv_capacity(small_scale, scaling_factors)],
            "y-large-scale-wind-gw": [wind_capacity(large_scale, scaling_factors)],
            "y-small-scale-wind-gw": [wind_capacity(small_scale, scaling_factors)],
            "y-large-scale-storage-gw": [storage_capacity_power(large_scale, scaling_factors)],
            "y-small-scale-storage-gw": [storage_capacity_power(small_scale, scaling_factors)],
            "y-large-scale-storage-gwh": [storage_capacity_energy(large_scale, scaling_factors)],
            "y-small-scale-storage-gwh": [storage_capacity_energy(small_scale, scaling_factors)],
            "y-large-scale-transmission-gwkm": [transmission_capacity(large_scale, scaling_factors)],
            "y-small-scale-transmission-gwkm": [transmission_capacity(small_scale, scaling_factors)],
            "y-large-scale-shed-load-gwh": [shed_load(large_scale, scaling_factors)],
            "y-small-scale-shed-load-gwh": [shed_load(small_scale, scaling_factors)]
        },
        index=[experiment_id]
    ).to_csv(path_to_output, sep="\t", index=True, header=True)


def system_cost(result, scaling_factors):
    return (result.results.cost
                  .squeeze("costs")
                  .sum("loc_techs_cost")
                  .item() / scaling_factors["monetary"])


def pv_capacity(result, scaling_factors):
    return (result.get_formatted_array("energy_cap")
                  .sel(techs=PV_TECHS)
                  .sum()
                  .item() / scaling_factors["power"] * ENERGY_SCALING_FACTOR)


def wind_capacity(result, scaling_factors):
    return (result.get_formatted_array("energy_cap")
                  .sel(techs=WIND_TECHS)
                  .sum()
                  .item() / scaling_factors["power"] * ENERGY_SCALING_FACTOR)


def storage_capacity_power(result, scaling_factors):
    return (result.results
                  .energy_cap
                  .sel(loc_techs=result.inputs.loc_techs_store)
                  .sum()
                  .item() / scaling_factors["power"] * ENERGY_SCALING_FACTOR)


def storage_capacity_energy(result, scaling_factors):
    return (result.results
                  .storage_cap
                  .sum()
                  .item() / scaling_factors["power"] * ENERGY_SCALING_FACTOR)


def transmission_capacity(result, scaling_factors):
    energy_cap = result.get_formatted_array("energy_cap") / scaling_factors["power"]
    trans_cap_mw = energy_cap.sel(techs=["ac_transmission" in str(tech) for tech in energy_cap.techs])
    if "distance" in result._model_data:
        distance_km = result.get_formatted_array("distance") / 1000
        return (trans_cap_mw * distance_km).sum().item() * ENERGY_SCALING_FACTOR / 2 # counting each line twice
    else:
        assert trans_cap_mw.sum() == 0
        return 0


def shed_load(result, scaling_factors):
    gen = result.get_formatted_array("carrier_prod") / scaling_factors["power"]
    if "load_shedding" in gen.techs:
        return gen.sel(techs="load_shedding").sum().item() * ENERGY_SCALING_FACTOR
    else:
        return 0


if __name__ == "__main__":
    determine_y(
        path_to_results_of_large_scale=snakemake.input.large_scale,
        path_to_results_of_small_scale=snakemake.input.small_scale,
        scaling_factors=snakemake.params.scaling_factors,
        experiment_id=snakemake.params.experiment_id,
        path_to_output=snakemake.output[0]
    )
