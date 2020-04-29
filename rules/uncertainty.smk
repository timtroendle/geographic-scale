"""Handle uncertainties of input parameters.

I am trying to investigate the impact of uncertain parameters on the difference
in total system cost of a large scale system and a small scale system.

In here, I assume the model to be a function y(x) that maps from parameter values x
to output value y. x is a vector, y is the difference in total system cost between
a large scale system and a small scale system.
"""
import pandas as pd

localrules: all_uncertainty, x, weather_diff, normal_diff, weather_diff_diff


rule all_uncertainty:
    message: "Perform entire uncertainty analysis."
    input:
        weather = "build/output/{resolution}/uncertainty/weather-diff-diff.txt".format(resolution=config["weather-uncertainty"]["resolution"]["space"]),


def ed_index(number_samples):
    """Creates complete index for experimental design based on number of samples."""
    return [f"_{x:03d}" for x in range(int(number_samples))]


rule experimental_design:
    message: "{wildcards.experiment}--{wildcards.number_samples}: Sample input for uncertainty analysis {wildcards.number_samples} times using UQLab."
    input:
        src = "src/uncertainty/experimental_design.m",
        parameters = "build/uncertainty/uncertain-parameters.csv"
    output: "build/uncertainty/{experiment}--{number_samples}/experimental-design.csv"
    run:
        import pandas as pd

        shell("""matlab -nodisplay -nodesktop -r \"cd src/uncertainty/; experimental_design('../../{input[1]}', '../../{output}', {wildcards.number_samples}); exit; \" """)

        ed = pd.read_csv(output[0])
        ed.index = ed_index(wildcards.number_samples)
        ed.to_csv(output[0], index=True, header=True)


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


rule biofuel_costs:
    message: "Merge biofuel cost scenarios."
    input:
        low = eurocalliope("build/data/{resolution}/biofuel/low/costs-eur-per-mwh.csv"),
        medium = eurocalliope("build/data/{resolution}/biofuel/medium/costs-eur-per-mwh.csv"),
        high = eurocalliope("build/data/{resolution}/biofuel/high/costs-eur-per-mwh.csv")
    params:
        efficiency = config["parameters"]["biofuel-efficiency"],
        om_prod = 3.1
    output: "build/uncertainty/{resolution}/c_bio.csv"
    run:
        import pandas as pd

        (
            pd
            .Series(
                index=["low", "medium", "high"],
                data=[pd.read_csv(x, header=None).iloc[0, 0] for x in [input.low, input.medium, input.high]],
                name="fuel_cost_eur_per_mwh"
            )
            .div(params.efficiency)
            .add(params.om_prod)
            .to_csv(output[0], header=True, index=True)
        )


rule x:
    message: "{wildcards.experiment}--{wildcards.number_samples}: Preprocess x {wildcards.id} to something usable by Calliope."
    input:
        src = "src/uncertainty/x.py",
        biofuel_potentials = rules.biofuel_availability.output[0],
        biofuel_costs = rules.biofuel_costs.output[0],
        experimental_design = rules.experimental_design.output[0]
    params:
        parameter_definitions = config["uncertainty"]["parameters"],
        scaling_factors = config["scaling-factors"]
    output: "build/uncertainty/{resolution}/{experiment}--{number_samples}/{id}-x.yaml"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/x.py"


rule national_uncertainty_run:
    message:
        "{wildcards.experiment}--{wildcards.number_samples}: Experiment {wildcards.id} using scenario {wildcards.scenario} and national resolution."
    input:
        model = "build/model/national/model.yaml",
        parameters = "build/uncertainty/national/{experiment}--{number_samples}/{id}-x.yaml"
    params:
        subset_time = config["uncertainty"]["subset_time"],
        time_resolution = config["uncertainty"]["resolution"]["time"],
    output: protected("build/output/national/runs/uncertainty/{experiment}--{number_samples}/{scenario}/{id}.nc")
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
        "{wildcards.experiment}--{wildcards.number_samples}: Experiment {wildcards.id} using scenario {wildcards.scenario} and regional resolution."
    input:
        model = "build/model/regional/model.yaml",
        parameters = "build/uncertainty/regional/{experiment}--{number_samples}/{id}-x.yaml"
    params:
        subset_time = config["uncertainty"]["subset_time"],
        time_resolution = config["uncertainty"]["resolution"]["time"],
    output: protected("build/output/regional/runs/uncertainty/{experiment}--{number_samples}/{scenario}/{id}.nc")
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
    message: "{wildcards.experiment}--{wildcards.number_samples}: Calculate y {wildcards.id} of {wildcards.scenario} on {wildcards.resolution} resolution."
    input:
        src = "src/uncertainty/y.py",
        result = "build/output/{resolution}/runs/uncertainty/{experiment}--{number_samples}/{scenario}/{id}.nc"
    output: "build/output/{resolution}/uncertainty/{experiment}--{number_samples}/{scenario}/{id}-y.tsv"
    params:
        scaling_factors = config["scaling-factors"],
        experiment_id = lambda wildcards: wildcards.id
    conda: "../envs/calliope.yaml"
    script: "../src/uncertainty/y.py"


def ys_in_experimental_design(wildcards):
    # returns paths of all ys of the experimental design
    path_to_ys = "build/output/{resolution}/uncertainty/{experiment}--{number_samples}/{scenario}/{{id}}-y.tsv".format(**wildcards)
    experimental_design_index = ed_index(wildcards["number_samples"])
    return [path_to_ys.format(id=idx) for idx in experimental_design_index]


rule xy:
    message: "{wildcards.experiment}--{wildcards.number_samples}: Gather all y of {wildcards.scenario} on {wildcards.resolution} resolution and combine with x."
    input:
        x = rules.experimental_design.output[0],
        y = ys_in_experimental_design
    output: "build/output/{resolution}/uncertainty/{experiment}--{number_samples}/{scenario}/xy.csv"
    run:
        import pandas as pd

        x = pd.read_csv(input.x, index_col=0)
        all_ys = pd.concat(
            [pd.read_csv(path_to_y, sep="\t", index_col=0) for path_to_y in input.y],
            axis="index"
        )
        all_data = pd.concat(
            [x, all_ys],
            axis="columns"
        ).to_csv(
            output[0],
            index=True,
            header=True
        )


rule all_experimental_designs:
    message: "Create all experiment designs."
    input:
        low_fidelity = "build/uncertainty/{experiment}--{number_samples}/experimental-design.csv".format(
            experiment=config["uncertainty"]["experiment-parameters"]["name"],
            number_samples=config["uncertainty"]["experiment-parameters"]["number-samples"]["low-fidelity"],
        ),
        high_fidelity = "build/uncertainty/{experiment}--{number_samples}/experimental-design.csv".format(
            experiment=config["uncertainty"]["experiment-parameters"]["name"],
            number_samples=config["uncertainty"]["experiment-parameters"]["number-samples"]["high-fidelity"],
        )


rule all_experiments:
    message: "Run all experiments and merge results."
    input:
        src = "src/uncertainty/merge_xys.py",
        mf_lf_1 = "build/output/{resolution}/uncertainty/{experiment}--{number_samples}/{scenario}/xy.csv".format(
            experiment=config["uncertainty"]["experiment-parameters"]["name"],
            number_samples=config["uncertainty"]["experiment-parameters"]["number-samples"]["low-fidelity"],
            scenario=config["uncertainty"]["experiment-parameters"]["multi-fidelity-1"]["scenario"],
            resolution=config["uncertainty"]["experiment-parameters"]["multi-fidelity-1"]["low-fidelity"]
        ),
        mf_hf_1 = "build/output/{resolution}/uncertainty/{experiment}--{number_samples}/{scenario}/xy.csv".format(
            experiment=config["uncertainty"]["experiment-parameters"]["name"],
            number_samples=config["uncertainty"]["experiment-parameters"]["number-samples"]["high-fidelity"],
            scenario=config["uncertainty"]["experiment-parameters"]["multi-fidelity-1"]["scenario"],
            resolution=config["uncertainty"]["experiment-parameters"]["multi-fidelity-1"]["high-fidelity"]
        ),
        mf_lf_2 = "build/output/{resolution}/uncertainty/{experiment}--{number_samples}/{scenario}/xy.csv".format(
            experiment=config["uncertainty"]["experiment-parameters"]["name"],
            number_samples=config["uncertainty"]["experiment-parameters"]["number-samples"]["low-fidelity"],
            scenario=config["uncertainty"]["experiment-parameters"]["multi-fidelity-2"]["scenario"],
            resolution=config["uncertainty"]["experiment-parameters"]["multi-fidelity-1"]["low-fidelity"]
        ),
        mf_hf_2 = "build/output/{resolution}/uncertainty/{experiment}--{number_samples}/{scenario}/xy.csv".format(
            experiment=config["uncertainty"]["experiment-parameters"]["name"],
            number_samples=config["uncertainty"]["experiment-parameters"]["number-samples"]["high-fidelity"],
            scenario=config["uncertainty"]["experiment-parameters"]["multi-fidelity-2"]["scenario"],
            resolution=config["uncertainty"]["experiment-parameters"]["multi-fidelity-2"]["high-fidelity"]
        ),
        sf = "build/output/{resolution}/uncertainty/{experiment}--{number_samples}/{scenario}/xy.csv".format(
            experiment=config["uncertainty"]["experiment-parameters"]["name"],
            number_samples=config["uncertainty"]["experiment-parameters"]["number-samples"]["low-fidelity"],
            scenario=config["uncertainty"]["experiment-parameters"]["single-fidelity"]["scenario"],
            resolution=config["uncertainty"]["experiment-parameters"]["single-fidelity"]["resolution"]
        )
    output:
        mf_lf = 'build/uncertainty/{}/xy-mf-lf.csv'.format(config["uncertainty"]["experiment-parameters"]["name"]),
        mf_hf = 'build/uncertainty/{}/xy-mf-hf.csv'.format(config["uncertainty"]["experiment-parameters"]["name"]),
        sf = 'build/uncertainty/{}/xy-sf.csv'.format(config["uncertainty"]["experiment-parameters"]["name"])
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/merge_xys.py"


rule uncertainty_analysis:
    message: "Determine Sobol' indices and sample y using UQLab."
    input:
        src = "src/uncertainty/uq.m",
        mf_lf = rules.all_experiments.output.mf_lf,
        mf_hf = rules.all_experiments.output.mf_hf,
        sf = rules.all_experiments.output.sf,
        parameters = "build/uncertainty/uncertain-parameters.csv"
    output: # ALL OUTPUTS HARDCODED IN SCRIPT!!
        total_sobol_sf = "build/uncertainty/{}/total-sobol-sf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        first_sobol_sf = "build/uncertainty/{}/first-sobol-sf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        total_minus_first_sobol_sf = "build/uncertainty/{}/total-minus-first-sobol-sf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        total_sobol_mf = "build/uncertainty/{}/total-sobol-mf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        first_sobol_mf = "build/uncertainty/{}/first-sobol-mf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        total_minus_first_sobol_mf = "build/uncertainty/{}/total-minus-first-sobol-mf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        samples_sf = "build/uncertainty/{}/pce-samples-sf.csv".format(config["uncertainty"]["experiment-parameters"]["name"]),
        samples_mf = "build/uncertainty/{}/pce-samples-mf.csv".format(config["uncertainty"]["experiment-parameters"]["name"])
    shell:
        """
        matlab -nodisplay -nodesktop -r \"cd src/uncertainty/; uq('../../{{input.mf_lf}}', '../../{{input.mf_hf}}', '../../{{input.sf}}', '../../{{input.parameters}}', '../../build/uncertainty/{}/'); exit; \"
        """.format(config["uncertainty"]["experiment-parameters"]["name"])



rule test_uncertainty_runs:
    message: "{wildcards.experiment}--{wildcards.number_samples}: Run tests for uncertainty runs of {wildcards.scenario} on {wildcards.resolution} resolution."
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        "tests/test_constraints.py",
        "tests/test_assumptions.py",
        results = ys_in_experimental_design,
        biofuel_potentials = eurocalliope("build/data/{{resolution}}/biofuel/{bioscenario}/potential-mwh-per-year.csv".format(
            bioscenario=config["parameters"]["jrc-biofuel"]["scenario"]
        )),
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    output: "build/output/{resolution}/uncertainty/{wildcards.experiment}--{wildcards.number_samples}/{scenario}/test-report.html"
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
    output:
        publish = "build/output/overview-uncertain-parameters.csv",
        uqlab = "build/uncertainty/uncertain-parameters.csv"
    conda: "../envs/default.yaml"
    script: "../src/uncertainty/overview_parameters.py"
