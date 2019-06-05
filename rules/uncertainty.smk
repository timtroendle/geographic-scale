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
localrules: x


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
        """ # FIXME should fail when sub-optimal, but currently I cannot make it optimal
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}\
        --no_fail_when_infeasible\
        --override_dict="{{import: [{input.parameters}], \
                           model.subset_time: {params.subset_time}, \
                           model.time.function: resample, \
                           model.time.function_options: {{'resolution': '{params.time_resolution}'}}}}"
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




