rule run_national:
    message: "Run the model for scenario {wildcards.scenario} on national resolution."
    input:
        model = "build/model/national/model.yaml"
    params: subset_time = config["subset_time"]
    output: "build/output/national/{scenario}/results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario} \
            --override_dict "{{model.subset_time: {params.subset_time}}}"
        """


rule run_regional: # this is a copy of run_national which is necessary to have different cluster configs
    message: "Run the model for scenario {wildcards.scenario} on regional resolution."
    input:
        model = "build/model/regional/model.yaml"
    params: subset_time = config["subset_time"]
    output: "build/output/regional/{scenario}/results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario} \
            --override_dict "{{model.subset_time: {params.subset_time}}}"
        """


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        results = expand(
            "build/output/{resolution}/{{scenario}}/results.nc".format(resolution=config["resolution"]),
            scenario=config["scenarios"]
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
        shapes = "data/{resolution}/units.geojson",
        continental_shape = "data/continental/units.geojson",
        national_shapes = "data/national/units.geojson",
        regional_shapes = "data/regional/units.geojson",
        continental_result = "build/output/{resolution}/continental-autarky-100-continental-grid/results.nc",
        national_result = "build/output/{resolution}/national-autarky-100-national-grid/results.nc",
        regional_result = "build/output/{resolution}/regional-autarky-100-regional-grid/results.nc"
    params:
        scaling_factor_cost = config["scaling-factors"]["power"] / config["scaling-factors"]["monetary"]
    output: "build/output/{resolution}/map.png"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/map.py"
