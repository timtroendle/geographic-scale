from pathlib import Path

import calliope
import pandas as pd

ENERGY_SCALING_FACTOR = 1 / 1e3 # from MW(h) to GW(h)
PV_TECHS = ["open_field_pv", "roof_mounted_pv"]
WIND_TECHS = ["wind_offshore", "wind_onshore_monopoly", "wind_onshore_competing"]
HYDRO_TECHS = ["hydro_run_of_river"]
BIOMASS_TECHS = ["biofuel"]
CURTAILABLE_RE_TECHS = PV_TECHS + WIND_TECHS + HYDRO_TECHS
ALL_STORAGE_TECHS = ["battery", "biofuel", "hydro_reservoir", "hydrogen", "pumped_hydro"]
DEPLOYABLE_STORAGE_TECHS = ["battery", "hydrogen"]
DEMAND_TECH = "demand_elec"


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
            "Biomass capacity [GW]": [biomass_capacity(result, scaling_factors) for result in scenario_results.values()],
            "Storage capacity [GW]": [storage_capacity_power(result, scaling_factors) for result in scenario_results.values()],
            "Storage capacity [GWh]": [storage_capacity_energy(result, scaling_factors) for result in scenario_results.values()],
            "Transmission capacity [GW km]": [transmission_capacity(result, scaling_factors) for result in scenario_results.values()],
            "Curtailment [%]": [relative_curtailment(result, scaling_factors) for result in scenario_results.values()],
            "Load shed [‰]": [shed_load(result, scaling_factors) for result in scenario_results.values()]
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
    return capacity(result, PV_TECHS, scaling_factors)


def wind_capacity(result, scaling_factors):
    return capacity(result, WIND_TECHS, scaling_factors)


def biomass_capacity(result, scaling_factors):
    return capacity(result, BIOMASS_TECHS, scaling_factors)


def storage_capacity_power(result, scaling_factors):
    return capacity(result, DEPLOYABLE_STORAGE_TECHS, scaling_factors)


def capacity(result, techs, scaling_factors):
    return (result.get_formatted_array("energy_cap")
                  .sel(techs=techs)
                  .sum()
                  .item() / scaling_factors["power"] * ENERGY_SCALING_FACTOR)


def storage_capacity_energy(result, scaling_factors):
    storage_cap = result.get_formatted_array("storage_cap")
    assert set(storage_cap.techs.values) == set(ALL_STORAGE_TECHS)
    return (storage_cap.sel(techs=DEPLOYABLE_STORAGE_TECHS)
                       .sum(["techs", "locs"])
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
    resource = data.get_formatted_array("resource").sel(techs=CURTAILABLE_RE_TECHS)
    capacity = data.get_formatted_array("energy_cap").sel(techs=CURTAILABLE_RE_TECHS) / scaling_factors["power"]
    potential = (resource * capacity).sum(dim=["timesteps", "techs", "locs"])
    generated = (data.get_formatted_array("carrier_prod")
                     .sel(techs=CURTAILABLE_RE_TECHS)
                     .sum(dim=["timesteps", "techs", "locs"])) / scaling_factors["power"]
    return ((potential - generated) / potential).to_pandas().electricity * 100


def shed_load(result, scaling_factors):
    gen = result.get_formatted_array("carrier_prod")
    if "load_shedding" in gen.techs:
        shed = gen.sel(techs="load_shedding").sum().item()
        demand = result.get_formatted_array("carrier_con").sel(techs=DEMAND_TECH).sum().item() * -1
        return (shed / demand * 1000)
    else:
        return 0


if __name__ == "__main__":
    main(
        paths_to_results=snakemake.input.results,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0]
    )
