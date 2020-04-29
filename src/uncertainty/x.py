from itertools import chain

import pandas as pd


def preprocess_x(experimental_design, idx, scaling_factors, parameter_definitions,
                 path_to_biofuel_potentials, path_to_fuel_costs, path_to_output):
    parameter_definitions = {
        parameter_definition["short-name"]: parameter_definition
        for long_name, parameter_definition in parameter_definitions.items()
    }
    x = pd.read_csv(experimental_design, index_col=0).loc[idx].to_dict()
    # special case treatment for ac_bio
    biofuel_potentials = (pd
                          .read_csv(path_to_biofuel_potentials, index_col=0)
                          .rename(index=lambda idx: idx.replace(".", "-")))
    biofuel_costs = pd.read_csv(path_to_fuel_costs, index_col=0).iloc[:, 0]
    ac_bio = _calliope_parameters_for_ac_bio(
        value=x["ac_bio"],
        potentials=biofuel_potentials * scaling_factors["power"],
        costs=biofuel_costs * scaling_factors["monetary"] / scaling_factors["power"]
    )
    del parameter_definitions["ac_bio"]
    del x["ac_bio"]
    # all other cases
    others = [
        _calliope_parameters(parameter_definitions[short_name], value, scaling_factors)
        for short_name, value in x.items()
    ]
    x_preprocessed = chain(*(others + [ac_bio]))

    with open(path_to_output, "w") as f_calliope:
        for line in x_preprocessed:
            f_calliope.write(f"{line}\n")


def _calliope_parameters(parameter_definition, value, scaling_factors):
    value = _scale_value(value, parameter_definition, scaling_factors)
    return [
        f"{calliope_param}: {value * scale}"
        for calliope_param, scale in parameter_definition["calliope-params"].items()
    ]


def _calliope_parameters_for_ac_bio(value, potentials, costs):
    assert value >= 0
    assert value <= 1
    if value <= 0.5:
        potential = potentials["low"] + value * 2 * (potentials["medium"] - potentials["low"])
        cost = costs["low"] + value * 2 * (costs["medium"] - costs["low"])
    if value > 0.5:
        potential = potentials["medium"] + (value - 0.5) * 2 * (potentials["high"] - potentials["medium"])
        cost = costs["medium"] + (value - 0.5) * 2 * (costs["high"] - costs["medium"])
    assert (potential >= potentials["low"]).all()
    assert (potential <= potentials["high"]).all()
    assert cost >= costs.loc["high"]
    assert cost <= costs.loc["low"]
    resources = [
        f"locations.{loc}.techs.biofuel.constraints.resource: {res / 8760}"
        for loc, res in potential.to_dict().items()
    ]
    storage = [
        f"locations.{loc}.techs.biofuel.constraints.storage_cap_equals: {res / 2}"
        for loc, res in potential.to_dict().items()
    ]
    om_con = [f"techs.biofuel.costs.monetary.om_con: {cost}"]
    return om_con + resources + storage


def _map_name(name, parameter_definitions):
    return parameter_definitions[name]["calliope-param"]


def _scale_value(value, parameter_definition, scaling_factors):
    assert value >= parameter_definition["min"], parameter_definition["short-name"]
    assert value <= parameter_definition["max"], parameter_definition["short-name"]
    scaling_factor = scaling_factors.get(parameter_definition["scaling-factor"], 1)
    inverse_scaling_factor = scaling_factors.get(parameter_definition["inverse-scaling-factor"], 1)
    return value * scaling_factor / inverse_scaling_factor


if __name__ == "__main__":
    preprocess_x(
        experimental_design=snakemake.input.experimental_design,
        idx=snakemake.wildcards.id,
        scaling_factors=snakemake.params.scaling_factors,
        parameter_definitions=snakemake.params.parameter_definitions,
        path_to_biofuel_potentials=snakemake.input.biofuel_potentials,
        path_to_fuel_costs=snakemake.input.biofuel_costs,
        path_to_output=snakemake.output[0]
    )
