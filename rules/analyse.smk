rule run_national:
    message: "Run the model for scenario {wildcards.scenario} on national resolution."
    input:
        model = "build/model/national/model.yaml"
    params:
        subset_time = config["subset_time"],
        time_resolution = config["resolution"]["time"]
    output: "build/output/national/{scenario}/results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario} \
            --override_dict "{{model.subset_time: {params.subset_time}, model.time.function: resample, \
                               model.time.function_options: {{'resolution': '{params.time_resolution}'}}}}"
        """


rule run_regional: # this is a copy of run_national which is necessary to have different cluster configs
    message: "Run the model for scenario {wildcards.scenario} on regional resolution."
    input:
        model = "build/model/regional/model.yaml"
    params:
        subset_time = config["subset_time"],
        time_resolution = config["resolution"]["time"]
    output: "build/output/regional/{scenario}/results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """ # because runs take a lot of computation time, don't fail when suboptimal
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario} \
            --no_fail_when_infeasible\
            --override_dict "{{model.subset_time: {params.subset_time}, model.time.function: resample, \
                               model.time.function_options: {{'resolution': '{params.time_resolution}'}}}}"
        """


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        "tests/test_constraints.py",
        results = expand(
            "build/output/{resolution}/{{scenario}}/results.nc".format(resolution=config["resolution"]["space"]),
            scenario=config["scenarios"]
        ),
        biofuel_potentials = eurocalliope("build/data/{resolution}/biofuel-potential-mwh-per-year.csv".format(
            resolution=config["resolution"]["space"])
        )
    output: "build/logs/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"


rule plot_scenario_space:
    message: "Plot scenario space and results."
    input:
        src = "src/analyse/scenarios.py",
        results = expand("build/output/{{resolution}}/{scenario}/results.nc", scenario=config["scenarios"])
    params:
        scaling_factor_cost = config["scaling-factors"]["power"] / config["scaling-factors"]["monetary"]
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
    params:
        scaling_factor_cost = config["scaling-factors"]["power"] / config["scaling-factors"]["monetary"]
    output: "build/output/{resolution}/map.png"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/map.py"


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


rule overview_scenario_results:
    message: "Create table of key outputs of scenarios."
    input:
        src = "src/analyse/scenario_overview.py",
        results = expand("build/output/{{resolution}}/{scenario}/results.nc", scenario=config["scenarios"])
    params: scaling_factors = config["scaling-factors"]
    output: "build/output/{resolution}/overview-scenario-results.csv"
    conda: "../envs/calliope.yaml"
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
