import xarray as xr


def weather_diff(path_to_large_scale, path_to_small_scale, path_to_output):
    costs_large_scale = read_total_system_costs(path_to_large_scale)
    costs_small_scale = read_total_system_costs(path_to_small_scale)
    with open(path_to_output, "w") as f_output:
        f_output.write(str(costs_small_scale / costs_large_scale))


def read_total_system_costs(path_to_results):
    return (xr.open_dataset(path_to_results)["total_levelised_cost"]
              .squeeze(["costs", "carriers"])
              .item())


if __name__ == "__main__":
    weather_diff(
        path_to_large_scale=snakemake.input.large_scale,
        path_to_small_scale=snakemake.input.small_scale,
        path_to_output=snakemake.output[0]
    )
