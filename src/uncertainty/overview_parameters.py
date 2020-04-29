import pandas as pd


def uncertainty_parameters_overview(parameters, path_to_uqlab, path_to_publish):
    pd.DataFrame(
        data={
            "Name": [p["descriptive-name"] for p in parameters.values()],
            "Description": [p["description"] for p in parameters.values()],
            "Min": [p["min"] for p in parameters.values()],
            "Max": [p["max"] for p in parameters.values()],
            "Unit": [p["unit"] for p in parameters.values()],
            "Source": [p["source"] for p in parameters.values()],
        }
    ).to_csv(path_to_publish, index=False, header=True)
    pd.DataFrame(
        data={
            "parameter-name": [p for p in parameters.keys()],
            "short-name": [p["short-name"] for p in parameters.values()],
            "min": [p["min"] for p in parameters.values()],
            "max": [p["max"] for p in parameters.values()],
            "unit": [p["unit"] for p in parameters.values()],
            "source": [p["source"] for p in parameters.values()],
        }
    ).to_csv(path_to_uqlab, index=False, header=True)


if __name__ == "__main__":
    uncertainty_parameters_overview(
        parameters=snakemake.params.parameters,
        path_to_uqlab=snakemake.output.uqlab,
        path_to_publish=snakemake.output.publish
    )
