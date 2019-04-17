rule run_national:
    message: "Run the model for scenario {wildcards.scenario} on national resolution."
    input:
        model = "build/model/national/model.yaml"
    output: "build/output/national/{scenario}/results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        "calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}"


rule run_regional: # this is a copy of run_national which is necessary to have different cluster configs
    message: "Run the model for scenario {wildcards.scenario} on regional resolution."
    input:
        model = "build/model/regional/model.yaml"
    output: "build/output/regional/{scenario}/results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        "calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}"


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        results = expand("build/output/national/{scenario}/results.nc", scenario=config["scenarios"])
    output: "build/logs/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"


rule plot_scenario_space:
    message: "Plot scenario space and results."
    input:
        src = "src/analyse/scenarios.py",
        results = expand("build/output/national/{scenario}/results.nc", scenario=config["scenarios"])
    output: "build/scenario-space.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/scenarios.py"


rule plot_map:
    message: "Create map of results."
    input:
        src = "src/analyse/map.py",
        continental_shape = "data/units-continental.geojson",
        national_shapes = "data/units-national.geojson",
        regional_shapes = "data/units-regional.geojson",
        continental_result = "build/output/national/continental-autarky-100-continental-grid/results.nc",
        national_result = "build/output/national/national-autarky-100-national-grid/results.nc"
    output: "build/map.png"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/map.py"
