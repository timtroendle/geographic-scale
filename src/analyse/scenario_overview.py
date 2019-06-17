from pathlib import Path

import calliope
import pandas as pd

ENERGY_SCALING_FACTOR = 1 / 1e3 # from MW(h) to GW(h)
PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
WIND_TECHS = ["wind_offshore", "wind_onshore_monopoly", "wind_onshore_competing"]
HYDRO_TECHS = ["hydro_run_of_river"]
RE_TECHS = PV_TECHS + WIND_TECHS + HYDRO_TECHS


def main(paths_to_results, scaling_factors, path_to_output):
    scenario_results = {
        Path(path_to_results).parent.name: calliope.read_netcdf(path_to_results)
        for path_to_results in paths_to_results
    }

    data = pd.DataFrame(
        index=scenario_results.keys(),
        data={
            "PV capacity [GW]": [pv_capacity(result, scaling_factors) for result in scenario_results.values()],
            "Wind capacity [GW]": [wind_capacity(result, scaling_factors) for result in scenario_results.values()],
            "Storage capacity [GW]": [storage_capacity_power(result, scaling_factors) for result in scenario_results.values()],
            "Storage capacity [GWh]": [storage_capacity_energy(result, scaling_factors) for result in scenario_results.values()],
            "Transmission capacity [GW km]": [transmission_capacity(result, scaling_factors) for result in scenario_results.values()],
            "Curtailment [%]": [relative_curtailment(result, scaling_factors) for result in scenario_results.values()],
            "Load shed [GWh]": [shed_load(result, scaling_factors) for result in scenario_results.values()]
        }
    )
    data.index.name = "Scenario"
    data.to_csv(
        path_to_output,
        index=True,
        header=True,
        float_format="%.0f"
    )


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


def relative_curtailment(data, scaling_factors):
    resource = data.get_formatted_array("resource").sel(techs=RE_TECHS)
    capacity = data.get_formatted_array("energy_cap").sel(techs=RE_TECHS) / scaling_factors["power"]
    potential = (resource * capacity).sum(dim=["timesteps", "techs", "locs"])
    generated = (data.get_formatted_array("carrier_prod")
                     .sel(techs=RE_TECHS)
                     .sum(dim=["timesteps", "techs", "locs"])) / scaling_factors["power"]
    return ((potential - generated) / potential).to_pandas().electricity * 100


def shed_load(result, scaling_factors):
    gen = result.get_formatted_array("carrier_prod") / scaling_factors["power"]
    if "load_shedding" in gen.techs:
        return gen.sel(techs="load_shedding").sum().item() * ENERGY_SCALING_FACTOR
    else:
        return 0


if __name__ == "__main__":
    main(
        paths_to_results=snakemake.input.results,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
