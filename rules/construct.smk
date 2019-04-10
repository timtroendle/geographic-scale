subworkflow eurocalliope:
    workdir: "./euro-calliope"
    snakefile: "./euro-calliope/Snakefile"
    configfile: "./config/default.yaml"

configfile: "./config/default.yaml"

localrules: copy_euro_calliope, model
ruleorder: model > copy_euro_calliope


rule copy_euro_calliope:
    message: "Copy file ./build/model/{wildcards.definition_file}.{wildcards.suffix} from euro-calliope."
    input: eurocalliope("build/model/{definition_file}.{suffix}"),
    output: "build/model/{definition_file}.{suffix}"
    shell: "cp {input} {output}"


rule model:
    message: "Build entire model."
    input:
        "build/model/interest-rate.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        "build/model/link-techs.yaml",
        "build/model/national/locations.yaml",
        "build/model/national/link-all-neighbours.yaml",
        expand(
            "build/model/{resolution}/electricity-demand.csv",
            resolution=["national"]
        ),
        expand(
            "build/model/{resolution}/capacityfactors-{technology}.csv",
            technology=["open-field-pv", "rooftop-pv", "wind-offshore", "wind-onshore"],
            resolution=["national"]
        ),
        definition = "src/template/model.yaml"
    output:
        model = "build/model/model.yaml"
    shell:
        "cp {input.definition} {output}"
