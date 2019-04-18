def preprocess_x(x, scaling_factors, parameter_definitions, path_to_output):
    x_preprocessed = {
        _map_name(name, parameter_definitions): _scale_value(name, value, parameter_definitions, scaling_factors)
        for name, value in x.to_dict().items()
    }

    with open(path_to_output, "w") as f_calliope:
        for name, value in x_preprocessed.items():
            f_calliope.write(f"{name}: {value}\n")


def _map_name(name, parameter_definitions):
    return parameter_definitions[name]["calliope-param"]


def _scale_value(name, value, parameter_definitions, scaling_factors):
    parameter_definition = parameter_definitions[name]
    scaling_factor = scaling_factors[parameter_definition["scaling-factor"]]
    inverse_scaling_factor = scaling_factors[parameter_definition["inverse-scaling-factor"]]
    return value * scaling_factor / inverse_scaling_factor


if __name__ == "__main__":
    preprocess_x(
        x=snakemake.params.x,
        scaling_factors=snakemake.params.scaling_factors,
        parameter_definitions=snakemake.params.parameter_definitions,
        path_to_output=snakemake.output[0]
    )
