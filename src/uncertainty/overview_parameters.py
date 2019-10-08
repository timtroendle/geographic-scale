import pandas as pd


def uncertainty_parameters_overview(parameters, path_to_output):
    pd.DataFrame(
        data={
            "Name": [p["descriptive-name"] for p in parameters.values()],
            "Description": [p["description"] for p in parameters.values()],
            "Min": [p["min"] for p in parameters.values()],
            "Max": [p["max"] for p in parameters.values()],
            "Unit": [p["unit"] for p in parameters.values()],
            "Source": [p["source"] for p in parameters.values()],
        }
    ).to_csv(path_to_output, index=False, header=True)


if __name__ == "__main__":
    uncertainty_parameters_overview(
        parameters=snakemake.params.parameters,
        path_to_output=snakemake.output[0]
    )
