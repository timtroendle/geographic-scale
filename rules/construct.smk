subworkflow eurocalliope:
    workdir: "./euro-calliope"
    snakefile: "./euro-calliope/Snakefile"
    configfile: "./config/default.yaml"

configfile: "./config/default.yaml"

localrules: copy_euro_calliope, load_shedding, model
ruleorder: model > import_restrictions > grid_size_restrictions > load_shedding > copy_euro_calliope

LOCATIONS_WITH_LOAD_SHEDDING = "data/{resolution}/locations-with-loadshedding.txt"


rule copy_euro_calliope:
    message: "Copy file ./build/model/{wildcards.definition_file}.{wildcards.suffix} from euro-calliope."
    input: eurocalliope("build/model/{definition_file}.{suffix}"),
    output: "build/model/{definition_file}.{suffix}"
    shell: "cp {input} {output}"


rule import_restrictions:
    message: "Create import restriction overrides for {wildcards.resolution} resolution."
    input:
        src = "src/construct/import_restrictions.py",
        units = eurocalliope("build/data/{resolution}/units.geojson")
    params:
        restrictions = [0, 15, 30]
    conda: "../envs/geo.yaml"
    output: "build/model/{resolution}/import-restrictions.yaml"
    script: "../src/construct/import_restrictions.py"


rule grid_size_restrictions:
    message: "Create grid size restriction overrides for {wildcards.resolution} resolution."
    input:
        src = "src/construct/grid_size_restrictions.py",
        links = "build/model/{resolution}/link-all-neighbours.yaml",
        units = eurocalliope("build/data/{resolution}/units.geojson")
    conda: "../envs/geo.yaml"
    output: "build/model/{resolution}/grid-size-restrictions.yaml"
    script: "../src/construct/grid_size_restrictions.py"


def load_shedding_locations(wildcards):
    locations = config["load-shedding"][wildcards["resolution"]]
    if locations:
        return locations
    else:
        return []


rule load_shedding:
    message: "Create load shedding override."
    input:
        src = "src/construct/load_shedding.py"
    params: locations = load_shedding_locations
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/load-shedding.yaml"
    script: "../src/construct/load_shedding.py"


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
        "build/model/{resolution}/load-shedding.yaml",
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
