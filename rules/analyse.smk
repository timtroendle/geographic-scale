localrules: overview_scenario_results

RESAMPLE_OVERRIDE_DICT = {
    "model.time.function": "resample",
    "model.time.function_options": {"resolution": config["resolution"]["time"]}
}
CALLIOPE_OVERRIDE_DICT = {
    **config["calliope-parameters"],
    **RESAMPLE_OVERRIDE_DICT
}


rule run_national:
    message: "Run the model for scenario {wildcards.scenario} on national resolution."
    input:
        model = "build/model/national/model.yaml"
    params: override_dict = CALLIOPE_OVERRIDE_DICT
    output: "build/output/national/{scenario}/results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario} \
            --override_dict="{params.override_dict}"
        """


rule run_regional: # this is a copy of run_national which is necessary to have different cluster configs
    message: "Run the model for scenario {wildcards.scenario} on regional resolution."
    input:
        model = "build/model/regional/model.yaml"
    params: override_dict = CALLIOPE_OVERRIDE_DICT
    output: "build/output/regional/{scenario}/results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """ # because runs take a lot of computation time, don't fail when suboptimal
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario} \
            --no_fail_when_infeasible\
            --override_dict="{params.override_dict}"
        """


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        "tests/test_constraints.py",
        "tests/test_assumptions.py",
        results = expand(
            "build/output/{{resolution}}/{scenario}/results.nc",
            scenario=config["scenarios"]
        ),
        biofuel_potentials = eurocalliope("build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        )),
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    output: "build/logs/{resolution}/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"


rule aggregated_results:
    message: "Create csv overview over all results."
    input:
        src = "src/analyse/aggregation.py",
        scenarios = expand("build/output/{{resolution}}/{scenario}/results.nc", scenario=config["scenarios"]),
        units = eurocalliope("build/data/{resolution}/units.geojson")
    params: scaling_factors = config["scaling-factors"]
    conda: "../envs/geo.yaml"
    output: "build/output/{resolution}/aggregation.csv"
    script: "../src/analyse/aggregation.py"


rule time_aggregated_results:
    message: "Create NetCDF overview over all results."
    input:
        src = "src/analyse/time_aggregation.py",
        scenarios = expand("build/output/{{resolution}}/{scenario}/results.nc", scenario=config["scenarios"]),
        units = eurocalliope("build/data/{resolution}/units.geojson")
    params: scaling_factors = config["scaling-factors"]
    conda: "../envs/geo.yaml"
    output: "build/output/{resolution}/aggregation.nc"
    script: "../src/analyse/time_aggregation.py"


rule plot_cost:
    message: "Plot overview over cost."
    input:
        src = "src/analyse/cost.py",
        results = rules.time_aggregated_results.output[0]
    output:
        base = "build/output/{resolution}/cost.{plot_suffix}",
        special = "build/output/{resolution}/cost-special-cases.{plot_suffix}",
    conda: "../envs/default.yaml"
    script: "../src/analyse/cost.py"


rule plot_map_cost:
    message: "Create map of cost."
    input:
        src = "src/analyse/map_cost.py",
        shapes = eurocalliope("build/data/{resolution}/units.geojson"),
        national_shapes = eurocalliope("build/data/national/units.geojson"),
        regional_shapes = eurocalliope("build/data/regional/units.geojson"),
        results = rules.time_aggregated_results.output[0]
    params:
        connected_regions = config["connected-regions"]
    output: "build/output/{resolution}/map-cost.png"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/map_cost.py"


rule plot_map_energy:
    message: "Create map of generation and transmission."
    input:
        src = "src/analyse/map_energy.py",
        shapes = eurocalliope("build/data/{resolution}/units.geojson"),
        results = rules.time_aggregated_results.output[0]
    params:
        connected_regions = config["connected-regions"],
        scenarios = ["continental-autarky-100-continental-grid",
                     "national-autarky-100-continental-grid",
                     "regional-autarky-100-continental-grid"]
    output: "build/output/{resolution}/map-energy.png"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/map_energy.py"


rule plot_network_map:
    message: "Create map of electricity grid."
    input:
        src = "src/analyse/network.py",
        units = eurocalliope("build/data/{resolution}/units.geojson"),
        results = "build/output/{resolution}/continental-autarky-100-continental-grid/results.nc"
    output: "build/output/{resolution}/network.{plot_suffix}"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/network.py"


rule plot_system_composition:
    message: "Create plot of system composition of all scenarios."
    input:
        src = "src/analyse/composition.py",
        results = rules.aggregated_results.output[0]
    output: "build/output/{resolution}/composition.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/composition.py"


rule plot_uncertainty_of_composition:
    message: "Create plot of uncertainty of system composition found in uncertainty analysis."
    input:
        src = "src/analyse/composition_uncertainty.py",
        large_scales = "data/pce-samples-continental-national-scales.csv",
        small_scale = "data/pce-samples-regional-scale.csv",
        aggregate = rules.time_aggregated_results.output[0]
    output: "build/output/{resolution}/composition-uncertainty.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/composition_uncertainty.py"


rule plot_uncertainty_of_cost:
    message: "Create plot of uncertainty of system cost found in uncertainty analysis."
    input:
        src = "src/analyse/cost_uncertainty.py",
        large_scales = "data/pce-samples-continental-national-scales.csv",
        small_scale = "data/pce-samples-regional-scale.csv",
        aggregate = rules.time_aggregated_results.output[0]
    output: "build/output/{resolution}/cost-uncertainty.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/cost_uncertainty.py"


rule plot_timeseries:
    message: "Create plot of timeseries on {wildcards.resolution} resolution."
    input:
        src = "src/analyse/timeseries.py",
        result = "build/output/{resolution}/regional-autarky-100-regional-grid/results.nc",
        units = eurocalliope("build/data/{resolution}/units.csv")
    params:
        connected_regions = config["connected-regions"],
        scaling_factors = config["scaling-factors"],
        resolution = config["report"]["timeseries-plot"]["resolution"]
    output:
        plot = "build/output/{resolution}/timeseries.{plot_suffix}",
        units = "build/output/{resolution}/timeseries-selection-{plot_suffix}.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/timeseries.py"


rule plot_sobol_indices:
    message: "Create heatmap of {wildcards.order} Sobol indices of all input/output combinations."
    input:
        src = "src/analyse/sobol.py",
        indices_cont_and_nat = "data/{order}-sobol-continental-national.csv",
        indices_reg = "data/{order}-sobol-regional.csv"
    params: parameters = config["uncertainty"]["parameters"]
    wildcard_constraints:
        order = "((total)|(first)|(total-minus-first))",
        extent = "((all)|(diff))"
    output: "build/output/{resolution}/{order}-sobol-{extent}.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/sobol.py"


rule overview_scenario_results:
    message: "Create table of key outputs of scenarios."
    input:
        src = "src/analyse/scenario_overview.py",
        results = rules.aggregated_results.output[0],
    output:
        table1 = "build/output/{resolution}/overview-scenario-results-1.csv",
        table2 = "build/output/{resolution}/overview-scenario-results-2.csv"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/scenario_overview.py"


rule overview_cost_assumptions:
    message: "Create table of key cost assumptions."
    input:
        src = "src/analyse/cost_assumptions.py",
        model = "build/output/{resolution}/regional-autarky-70-continental-grid/results.nc"
    params: scaling_factors = config["scaling-factors"]
    output: "build/output/{resolution}/overview-cost-assumptions.csv"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/cost_assumptions.py"
