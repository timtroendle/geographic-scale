"""Rules to assess where autarky is impossible, based on feasibility tests.

locations, a single location model will be infeasible. The rules in this file can run
infeasibility tests for all locations.

Because the test is computationally expensive, these rules are not in-the-loop. Instead,
they are run manually and the result is then checked in into `config/default.yaml`.
"""
localrules: single_location_model
ruleorder: single_location_model > copy_euro_calliope


def locations(resolution):
    with open(f"euro-calliope/build/model/{resolution}/electricity-demand.csv") as f_elec:
        columns = f_elec.readline()
    return columns.strip().split(",")[1:]


rule single_location_model:
    message: "Build model for single location runs on resolution {wildcards.resolution}."
    input:
        "build/model/interest-rate.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        "build/model/link-techs.yaml",
        "build/model/{resolution}/locations.yaml",
        "build/model/{resolution}/electricity-demand.csv",
        expand(
            "build/model/{{resolution}}/capacityfactors-{technology}.csv",
            technology=["open-field-pv", "rooftop-pv", "wind-offshore", "wind-onshore",
                        "hydro-ror", "hydro-reservoir-inflow"],
        ),
        definition = "src/template/model-single-location.yaml"
    output: "build/model/{resolution}/model-single-location.yaml"
    shell: "cp {input.definition} {output}"


rule run_single_location:
    message: "Run the model for single location {wildcards.location}."
    input:
        src = "src/impossible/run_single_location.py",
        model = "build/model/{resolution}/model-single-location.yaml"
    params:
        subset_time = config["subset_time"],
        time_resolution = "1H", #config["resolution"]["time"],
        all_locations = lambda wildcards: locations(wildcards["resolution"])
    output: "build/output/{resolution}/single-location/{location}.txt"
    log: "build/logs/single-location/{resolution}/{location}"
    conda: "../envs/calliope.yaml"
    script: "../src/impossible/run_single_location.py"
