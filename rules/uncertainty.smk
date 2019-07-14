"""Handle uncertainties of input parameters.

I am trying to investigate the impact of uncertain parameters on the difference
in total system cost of a large scale system and a small scale system.

In here, I assume the model to be a function y(x) that maps from parameter values x
to output value y. x is a vector, y is the difference in total system cost between
a large scale system and a small scale system.
"""
import pandas as pd

experiment_design = pd.read_csv(config["uncertainty"]["experiment-design"], index_col=0, sep="\t")
experiment_design.index = experiment_design.index.astype(str)
localrules: x, weather_diff, normal_diff, weather_diff_diff


rule x:
    message: "Preprocess x {wildcards.id} to something usable by Calliope."
    input: src = "src/uncertainty/x.py"
    params:
        x = lambda wildcards: experiment_design.loc[str(wildcards.id), :],
        parameter_definitions = config["uncertainty"]["parameters"],
        scaling_factors = config["scaling-factors"]
    output: "build/uncertainty/{id}-x.yaml"
    script: "../src/uncertainty/x.py"


rule uncertainty_run:
    message:
        "Experiment {{wildcards.id}} using scenario {{wildcards.scenario}} and {} resolution.".format(
            config["uncertainty"]["resolution"]["space"]
        )
    input:
        model = "build/model/{}/model.yaml".format(config["uncertainty"]["resolution"]["space"]),
        parameters = rules.x.output
    params:
        subset_time = config["uncertainty"]["subset_time"],
        time_resolution = config["uncertainty"]["resolution"]["time"],
    output: "build/uncertainty/{id}--{scenario}-results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}\
        --override_dict="{{import: [{input.parameters}], \
                           model.subset_time: {params.subset_time}, \
                           model.time.function: resample, \
                           model.time.function_options: {{'resolution': '{params.time_resolution}'}}}}"
        """


rule weather_run:
    message:
        "Run weather sensitivity using scenario {wildcards.scenario} and {wildcards.resolution} resolution."
    input:
        model = "build/model/{resolution}/model.yaml",
        elec_ts = "build/model/{resolution}/uncertainty/electricity-demand.csv",
        ror_ts = "build/model/{resolution}/uncertainty/capacityfactors-hydro-ror.csv",
        reservoir_ts = "build/model/{resolution}/uncertainty/capacityfactors-hydro-reservoir-inflow.csv"
    params:
        start_year = config["weather-uncertainty"]["start-year"],
        final_year = config["year"],
        time_resolution = config["weather-uncertainty"]["resolution"]["time"],
    output: "build/output/{resolution}/uncertainty/weather-{scenario}.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}\
        --override_dict="{{model.subset_time: [{params.start_year}-01-01, {params.final_year}-12-31], \
                           model.time.function: resample, \
                           model.time.function_options: {{'resolution': '{params.time_resolution}'}}, \
                           techs.demand_elec.constraints.resource: file=./uncertainty/electricity-demand.csv, \
                           techs.hydro_run_of_river.constraints.resource: file=./uncertainty/capacityfactors-hydro-ror.csv, \
                           techs.hydro_reservoir.constraints.resource: file=./uncertainty/capacityfactors-hydro-reservoir-inflow.csv, \
                           }}"
        """


rule y:
    message: "Calculate y for experiment {wildcards.id}."
    input:
        src = "src/uncertainty/y.py",
        large_scale = "build/uncertainty/{{id}}--{}-results.nc".format(config["uncertainty"]["scenarios"]["large-scale"]),
        small_scale = "build/uncertainty/{{id}}--{}-results.nc".format(config["uncertainty"]["scenarios"]["small-scale"])
    output: "build/uncertainty/{id}-y.tsv"
    params:
        scaling_factors = config["scaling-factors"],
        experiment_id = lambda wildcards: wildcards.id
    conda: "../envs/calliope.yaml"
    script: "../src/uncertainty/y.py"


rule xy:
    message: "Gather all y and combine with x."
    input:
        y = expand("build/uncertainty/{id}-y.tsv", id=experiment_design.index)
    output: "build/uncertainty/xy.csv"
    params: x = experiment_design
    run:
        import pandas as pd

        all_ys = pd.concat(
            [pd.read_csv(path_to_y, sep="\t", index_col=0) for path_to_y in input.y],
            axis="index"
        )
        all_data = pd.concat(
            [params.x, all_ys],
            axis="columns"
        ).to_csv(
            output[0],
            index=True,
            header=True
        )


rule repeat_timeseries:
    message: "Take single year timeseries {wildcards.timeseries} and repeat it from {params.start} onwards."
    input:
        src = "src/uncertainty/repeat.py",
        timeseries = "build/model/{resolution}/{timeseries}.csv"
    params: start = config["weather-uncertainty"]["start-year"]
    output: "build/model/{resolution}/uncertainty/{timeseries}.csv"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/repeat.py"


rule weather_diff:
    message: "Determine cost difference between large scale and small scale system for long weather."
    input:
        src = "src/uncertainty/weather_diff.py",
        large_scale = "build/output/{{resolution}}/uncertainty/weather-{}.nc".format(config["weather-uncertainty"]["scenarios"]["large-scale"]),
        small_scale = "build/output/{{resolution}}/uncertainty/weather-{}.nc".format(config["weather-uncertainty"]["scenarios"]["small-scale"])
    output: "build/output/uncertainty/{resolution}/weather-diff.txt"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/weather_diff.py"


rule normal_diff:
    message: "Determine cost difference between large scale and small scale system for resolution {wildcards.resolution}."
    input:
        src = "src/uncertainty/weather_diff.py",
        large_scale = "build/output/{{resolution}}/{}/results.nc".format(config["weather-uncertainty"]["scenarios"]["large-scale"]),
        small_scale = "build/output/{{resolution}}/{}/results.nc".format(config["weather-uncertainty"]["scenarios"]["small-scale"])
    output: "build/output/uncertainty/{resolution}/normal-diff.txt"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/weather_diff.py"


rule weather_diff_diff:
    message: "Determine the diff of weather runs to normal runs."
    input:
        normal = rules.normal_diff.output[0],
        weather = rules.weather_diff.output[0]
    output: "build/output/uncertainty/{resolution}/weather-diff-diff.txt"
    run:
        with open(input.normal, "r") as f_normal:
            normal_diff = float(f_normal.readline())
        with open(input.weather, "r") as f_weather:
            weather_diff = float(f_weather.readline())
        diff_diff = abs(normal_diff - weather_diff) / normal_diff
        with open(output[0], "w") as f_output:
            f_output.write(str(diff_diff))
