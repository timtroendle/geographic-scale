rule run:
    message: "Run the model for scenario {wildcards.scenario}."
    input:
        model = rules.model.output.model
    output: "build/output/{scenario}/results.nc"
    conda: "../envs/calliope.yaml"
    shell:
        "calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}"


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        results = expand("build/output/{scenario}/results.nc", scenario=config["scenarios"])
    output: "build/logs/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"


rule plot_scenario_space:
    message: "Plot scenario space and results."
    input:
        src = "src/analyse/scenarios.py",
        results = expand("build/output/{scenario}/results.nc", scenario=config["scenarios"])
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
        continental_result = "build/output/continental-autarky-100-continental-grid/results.nc",
        national_result = "build/output/national-autarky-100-national-grid/results.nc"
    output: "build/map.png"
    conda: "../envs/geo.yaml"
    script: "../src/analyse/map.py"
