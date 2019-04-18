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


rule x:
    message: "Preprocess x {wildcards.id} to something usable by Calliope."
    input: src = "src/uncertainty/x.py"
    params:
        x = lambda wildcards: experiment_design.loc[str(wildcards.id), :],
        parameter_definitions = config["uncertainty"]["parameters"],
        scaling_factors = config["scaling-factors"]
    output: "build/uncertainty/{id}-x.yaml"
    script: "../src/uncertainty/x.py"


rule model_run:
    message:
        "Experiment {{wildcards.id}} using scenario {{wildcards.scenario}} and {} resolution.".format(
            config["uncertainty"]["resolution"]
        )
    input:
        model = "build/model/{}/model.yaml".format(config["uncertainty"]["resolution"]),
        parameters = rules.x.output
    output: "build/uncertainty/{id}--{scenario}-results.nc"
    params: resolution = config["uncertainty"]["resolution"]
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}\
        --override_dict="{{import: [{input.parameters}] }}"
        """


rule y:
    message: "Calculate y for experiment {wildcards.id}."
    input:
        src = "src/uncertainty/y.py",
        large_scale = "build/uncertainty/{id}--continental-autarky-100-continental-grid-results.nc",
        small_scale = "build/uncertainty/{id}--regional-autarky-100-continental-grid-results.nc"
    output: "build/uncertainty/{id}-y.tsv"
    params:
        scaling_factors = config["scaling-factors"],
        experiment_id = lambda wildcards: wildcards.id
    conda: "../envs/default.yaml"
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




