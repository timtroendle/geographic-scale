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
localrules: all_uncertainty, x, weather_diff, normal_diff, weather_diff_diff


rule all_uncertainty:
    message: "Perform entire uncertainty analysis."
    input:
        weather = "build/output/{resolution}/uncertainty/weather-diff-diff.txt".format(resolution=config["weather-uncertainty"]["resolution"]["space"]),


rule biofuel_availability:
    message: "Merge biofuel availability scenarios."
    input:
        low = eurocalliope("build/data/{resolution}/biofuel/low/potential-mwh-per-year.csv"),
        medium = eurocalliope("build/data/{resolution}/biofuel/medium/potential-mwh-per-year.csv"),
        high = eurocalliope("build/data/{resolution}/biofuel/high/potential-mwh-per-year.csv")
    params: efficiency = config["parameters"]["biofuel-efficiency"]
    output: "build/uncertainty/{resolution}/a_bio.csv"
    run:
        import pandas as pd

        pd.DataFrame({
            "low": pd.read_csv(input.low, index_col=0).iloc[:, 0],
            "medium": pd.read_csv(input.medium, index_col=0).iloc[:, 0],
            "high": pd.read_csv(input.high, index_col=0).iloc[:, 0],
        }).mul(params.efficiency).to_csv(output[0], header=True, index=True)


rule x:
    message: "Preprocess x {wildcards.id} to something usable by Calliope."
    input:
        src = "src/uncertainty/x.py",
        biofuel = rules.biofuel_availability.output[0]
    params:
        x = lambda wildcards: experiment_design.loc[str(wildcards.id), :],
        parameter_definitions = config["uncertainty"]["parameters"],
        scaling_factors = config["scaling-factors"]
    output: "build/uncertainty/{resolution}/{id}-x.yaml"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/x.py"


rule national_uncertainty_run:
    message:
        "Experiment {wildcards.id} using scenario {wildcards.scenario} and national resolution."
    input:
        model = "build/model/national/model.yaml",
        parameters = "build/uncertainty/national/{id}-x.yaml"
    params:
        subset_time = config["uncertainty"]["subset_time"],
        time_resolution = config["uncertainty"]["resolution"]["time"],
    output: "build/output/national/uncertainty/{scenario}/{id}.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}\
        --override_dict="{{import: [{input.parameters}], \
                           model.subset_time: {params.subset_time}, \
                           model.time.function: resample, \
                           model.time.function_options: {{'resolution': '{params.time_resolution}'}}}}"
        """


rule regional_uncertainty_run:
    message:
        "Experiment {wildcards.id} using scenario {wildcards.scenario} and regional resolution."
    input:
        model = "build/model/regional/model.yaml",
        parameters = "build/uncertainty/regional/{id}-x.yaml"
    params:
        subset_time = config["uncertainty"]["subset_time"],
        time_resolution = config["uncertainty"]["resolution"]["time"],
    output: "build/output/regional/uncertainty/{scenario}/{id}.nc"
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
        elec_ts = "build/model/{resolution}/weather-uncertainty/electricity-demand.csv",
        ror_ts = "build/model/{resolution}/weather-uncertainty/capacityfactors-hydro-ror.csv",
        reservoir_ts = "build/model/{resolution}/weather-uncertainty/capacityfactors-hydro-reservoir-inflow.csv"
    params:
        override_dict = {
            **config["weather-uncertainty"]["calliope-parameters"],
            **{
                "model.subset_time": [f"{config['weather-uncertainty']['start-year']}-01-01", f"{config['year']}-12-31"],
                "model.time.function": "resample",
                "model.time.function_options": {'resolution': f'{config["weather-uncertainty"]["resolution"]["time"]}'},
                "techs.demand_elec.constraints.resource": "file=./weather-uncertainty/electricity-demand.csv",
                "techs.hydro_run_of_river.constraints.resource": "file=./weather-uncertainty/capacityfactors-hydro-ror.csv",
                "techs.hydro_reservoir.constraints.resource": "file=./weather-uncertainty/capacityfactors-hydro-reservoir-inflow.csv"
            }
        }
    output: "build/output/{resolution}/runs/uncertainty/weather-{scenario}.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}\
        --override_dict="{params.override_dict}"
        """


rule y:
    message: "Calculate y for experiment {wildcards.id} of {wildcards.scenario} on {wildcards.resolution} resolution."
    input:
        src = "src/uncertainty/y.py",
        result = "build/output/{resolution}/uncertainty/{scenario}/{id}.nc"
    output: "build/output/{resolution}/uncertainty/{scenario}/{id}-y.tsv"
    params:
        scaling_factors = config["scaling-factors"],
        experiment_id = lambda wildcards: wildcards.id
    conda: "../envs/calliope.yaml"
    script: "../src/uncertainty/y.py"


rule xy:
    message: "Gather all y of {wildcards.scenario} on {wildcards.resolution} resolution and combine with x."
    input:
        y = expand("build/output/{{resolution}}/uncertainty/{{scenario}}/{id}-y.tsv", id=experiment_design.index)
    output: "build/output/{resolution}/uncertainty/{scenario}/xy.csv"
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


rule test_uncertainty_runs:
    message: "Run tests for uncertainty runs of {wildcards.scenario} on {wildcards.resolution} resolution."
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        "tests/test_constraints.py",
        "tests/test_assumptions.py",
        results = expand(
            "build/output/{{resolution}}/uncertainty/{{scenario}}/{id}.nc",
            id=experiment_design.index
        ),
        biofuel_potentials = eurocalliope("build/data/{{resolution}}/biofuel/{bioscenario}/potential-mwh-per-year.csv".format(
            bioscenario=config["parameters"]["jrc-biofuel"]["scenario"]
        )),
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    output: "build/output/{resolution}/uncertainty/{scenario}/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"


rule repeat_timeseries:
    message: "Take single year timeseries {wildcards.timeseries} and repeat it from {params.start} onwards."
    input:
        src = "src/uncertainty/repeat.py",
        timeseries = "build/model/{resolution}/{timeseries}.csv"
    params: start = config["weather-uncertainty"]["start-year"]
    output: "build/model/{resolution}/weather-uncertainty/{timeseries}.csv"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/repeat.py"


rule weather_diff:
    message: "Determine cost difference between large scale and small scale system for long weather."
    input:
        src = "src/uncertainty/weather_diff.py",
        large_scale = "build/output/{{resolution}}/runs/uncertainty/weather-{}.nc".format(config["weather-uncertainty"]["scenarios"]["large-scale"]),
        small_scale = "build/output/{{resolution}}/runs/uncertainty/weather-{}.nc".format(config["weather-uncertainty"]["scenarios"]["small-scale"])
    output: "build/output/{resolution}/uncertainty/weather-diff.txt"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/weather_diff.py"


rule normal_diff:
    message: "Determine cost difference between large scale and small scale system for resolution {wildcards.resolution}."
    input:
        src = "src/uncertainty/weather_diff.py",
        large_scale = "build/output/{{resolution}}/runs/{}.nc".format(config["weather-uncertainty"]["scenarios"]["large-scale"]),
        small_scale = "build/output/{{resolution}}/runs/{}.nc".format(config["weather-uncertainty"]["scenarios"]["small-scale"])
    output: "build/output/{resolution}/uncertainty/normal-diff.txt"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/weather_diff.py"


rule weather_diff_diff:
    message: "Determine the diff of weather runs to normal runs."
    input:
        normal = rules.normal_diff.output[0],
        weather = rules.weather_diff.output[0]
    output: "build/output/{resolution}/uncertainty/weather-diff-diff.txt"
    run:
        with open(input.normal, "r") as f_normal:
            normal_diff = float(f_normal.readline())
        with open(input.weather, "r") as f_weather:
            weather_diff = float(f_weather.readline())
        diff_diff = abs(normal_diff - weather_diff) / normal_diff
        with open(output[0], "w") as f_output:
            f_output.write(str(diff_diff))


rule overview_uncertainty_parameters:
    message: "Create table of uncertainty parameters."
    input: src = "src/uncertainty/overview_parameters.py"
    params: parameters = config["uncertainty"]["parameters"]
    output: "build/output/{resolution}/overview-uncertain-parameters.csv"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/overview_parameters.py"
