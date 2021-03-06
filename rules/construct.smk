subworkflow eurocalliope:
    workdir: "./euro-calliope"
    snakefile: "./euro-calliope/Snakefile"
    configfile: "./config/default.yaml"

localrules: copy_euro_calliope, copy_resolution_specific_euro_calliope, model
ruleorder: model > import_restrictions > grid_size_restrictions > copy_euro_calliope > copy_resolution_specific_euro_calliope
wildcard_constraints:
    definition_file = "[^\/]*" # must not travers into directories


rule copy_euro_calliope:
    message: "Copy file {input[0]} from euro-calliope."
    input: eurocalliope("build/model/{definition_file}.{suffix}"),
    output: "build/model/{definition_file}.{suffix}"
    shell: "cp {input} {output}"


rule copy_resolution_specific_euro_calliope:
    message: "Copy file {input[0]} from euro-calliope."
    input: eurocalliope("build/model/{resolution}/{definition_file}.{suffix}"),
    output: "build/model/{resolution}/{definition_file}.{suffix}"
    shell: "cp {input} {output}"


rule import_restrictions:
    message: "Create import restriction overrides for {wildcards.resolution} resolution."
    input:
        src = "src/construct/import_restrictions.py",
        units = eurocalliope("build/data/{resolution}/units.geojson")
    params:
        restrictions = [0, 15, 30],
        connected_regions = config["connected-regions"]
    conda: "../envs/geo.yaml"
    output: "build/model/{resolution}/import-restrictions.yaml"
    script: "../src/construct/import_restrictions.py"


rule grid_size_restrictions:
    message: "Create grid size restriction overrides for {wildcards.resolution} resolution."
    input:
        src = "src/construct/grid_size_restrictions.py",
        links = "build/model/{resolution}/link-all-neighbours.yaml",
        units = eurocalliope("build/data/{resolution}/units.geojson")
    params:
        connected_regions = config["connected-regions"]
    conda: "../envs/geo.yaml"
    output: "build/model/{resolution}/grid-size-restrictions.yaml"
    script: "../src/construct/grid_size_restrictions.py"


rule model:
    message: "Build entire model on resolution {wildcards.resolution}."
    input:
        "build/model/interest-rate.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        "build/model/link-techs.yaml",
        "build/model/{resolution}/locations.yaml",
        "build/model/{resolution}/link-all-neighbours.yaml",
        "build/model/{resolution}/import-restrictions.yaml",
        "build/model/{resolution}/grid-size-restrictions.yaml",
        "build/model/{resolution}/electricity-demand.csv",
        expand(
            "build/model/{{resolution}}/capacityfactors-{technology}.csv",
            technology=["open-field-pv", "rooftop-pv", "wind-offshore", "wind-onshore",
                        "hydro-ror", "hydro-reservoir-inflow"],
        ),
        definition = "src/template/model.yaml"
    output:
        model = "build/model/{resolution}/model.yaml"
    shell:
        "cp {input.definition} {output}"
