import pandas as pd
import xarray as xr


def determine_y(path_to_results_of_large_scale, path_to_results_of_small_scale,
                scaling_factors, experiment_id, path_to_output):
    costs_large_scale = _excavate_cost_total_system(path_to_results_of_large_scale, scaling_factors)
    costs_small_scale = _excavate_cost_total_system(path_to_results_of_small_scale, scaling_factors)
    cost_diff = costs_small_scale - costs_large_scale
    pd.DataFrame(
        columns=["y_costs_large_scale", "y_costs_small_scale", "y_cost_diff"],
        data={
            "y_costs_large_scale": [costs_large_scale],
            "y_costs_small_scale": [costs_small_scale],
            "y_cost_diff": [cost_diff]
        },
        index=[experiment_id]
    ).to_csv(path_to_output, sep="\t", index=True, header=True)


def _excavate_cost_total_system(path_to_results, scaling_factors):
    return (xr.open_dataset(path_to_results)["cost"]
              .squeeze("costs")
              .sum("loc_techs_cost")
              .item() / scaling_factors["monetary"])


if __name__ == "__main__":
    determine_y(
        path_to_results_of_large_scale=snakemake.input.large_scale,
        path_to_results_of_small_scale=snakemake.input.small_scale,
        scaling_factors=snakemake.params.scaling_factors,
        experiment_id=snakemake.params.experiment_id,
        path_to_output=snakemake.output[0]
    )
