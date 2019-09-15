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


rule plot_scenario_space:
    message: "Plot scenario space and results."
    input:
        src = "src/analyse/scenarios.py",
        results = expand("build/output/{{resolution}}/{scenario}/results.nc", scenario=config["scenarios"])
    output: "build/output/{resolution}/scenario-space.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/scenarios.py"


rule plot_map:
    message: "Create map of results."
    input:
        src = "src/analyse/map.py",
        shapes = eurocalliope("build/data/{resolution}/units.geojson"),
        continental_shape = eurocalliope("build/data/continental/units.geojson"),
        national_shapes = eurocalliope("build/data/national/units.geojson"),
        regional_shapes = eurocalliope("build/data/regional/units.geojson"),
        continental_result = "build/output/{resolution}/continental-autarky-100-continental-grid/results.nc",
        national_result = "build/output/{resolution}/national-autarky-100-national-grid/results.nc",
        regional_result = "build/output/{resolution}/regional-autarky-100-regional-grid/results.nc"
    output: "build/output/{resolution}/map.png"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/map.py"


rule plot_network_map:
    message: "Create map of electricity grid."
    input:
        src = "src/analyse/network.py",
        units = eurocalliope("build/data/{resolution}/units.geojson"),
        results = "build/output/{resolution}/continental-autarky-100-continental-grid/results.nc"
    output: "build/output/{resolution}/network.png"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/network.py"


rule plot_electricity_flows:
    message: "Create map of electricity flows."
    input:
        src = "src/analyse/flows.py",
        units = eurocalliope("build/data/{resolution}/units.geojson"),
        countries = eurocalliope("build/data/national/units.geojson"),
        results = "build/output/{resolution}/continental-autarky-100-continental-grid/results.nc"
    params:
        scaling_factor = config["scaling-factors"]["power"]
    output: "build/output/{resolution}/flows.png"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/flows.py"


rule plot_system_composition:
    message: "Create plot of system composition of all scenarios."
    input:
        src = "src/analyse/composition.py",
        results = rules.aggregated_results.output[0]
    output: "build/output/{resolution}/composition.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/composition.py"


rule plot_variability_of_composition:
    message: "Create plot of variability of system composition found in uncertainty analysis."
    input:
        src = "src/analyse/composition_variability.py",
        xy = rules.xy.output[0], # FIXME should use data from surrogate model
        aggregate = rules.time_aggregated_results.output[0]
    output: "build/output/{resolution}/uncertainty/variability.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/composition_variability.py"


rule plot_timeseries:
    message: "Create plot of timeseries on regional resolution."
    input:
        src = "src/analyse/timeseries.py",
        result = "build/output/{resolution}/regional-autarky-100-regional-grid/results.nc",
        units = eurocalliope("build/data/{resolution}/units.csv")
    params:
        connected_regions = config["connected-regions"],
        scaling_factors = config["scaling-factors"],
        unit_lcoe = config["report"]["timeseries-plot"]["thresholds"]["unit-lcoe"],
        biofuel_lcoe = config["report"]["timeseries-plot"]["thresholds"]["biofuel-lcoe"],
        hydrogen_lcos = config["report"]["timeseries-plot"]["thresholds"]["hydrogen-lcos"],
        resolution = config["report"]["timeseries-plot"]["resolution"]
    output:
        plot = "build/output/{resolution}/timeseries.png",
        units = "build/output/{resolution}/timeseries-selection.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/timeseries.py"


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
