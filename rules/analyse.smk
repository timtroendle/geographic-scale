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
