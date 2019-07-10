import calliope


def run_single_location(path_to_yaml, path_to_log, location, all_locations,
                        subset_time, time_resolution, path_to_result):
    subset_time = subset_time[1:-1].split(",")
    subset_time = [x.strip()[1:-1] for x in subset_time]
    override_dict = {
        "model.subset_time": [subset_time[0], subset_time[1]],
        "model.time.function": "resample",
        "model.time.function_options": {"resolution": time_resolution},
        f"locations.{location}.exists": "True"
    }
    no_locations = {
        f"locations.{other_location}.exists": False
        for other_location in all_locations
        if other_location != location
    }
    model = calliope.Model(
        path_to_yaml,
        override_dict={**override_dict, **no_locations}
    )
    model.run_config['save_logs'] = path_to_log
    model.run()
    with open(path_to_result, "w") as f_result:
        f_result.write(model.results.termination_condition)


if __name__ == "__main__":
    run_single_location(
        path_to_yaml=snakemake.input.model,
        path_to_log=snakemake.log[0],
        location=snakemake.wildcards.location,
        all_locations=snakemake.params.all_locations,
        subset_time=snakemake.params.subset_time,
        time_resolution=snakemake.params.time_resolution,
        path_to_result=snakemake.output[0]
    )
