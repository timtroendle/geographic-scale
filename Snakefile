PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

configfile: "./config/default.yaml"
wildcard_constraints:
    resolution = "((continental)|(national)|(regional))", # supported spatial resolutions
    scenario = "({})".format("|".join([f"({scenario})" for scenario in config["scenarios"]])),
    plot_suffix = "((png)|(svg)|(tiff))"
onstart:
    shell("mkdir -p build/logs/uncertainty")
onsuccess:
    if "ifttt_apikey" in config.keys():
        trigger_ifttt(event_name="snakemake_succeeded", apikey=config["ifttt_apikey"])
onerror:
    if "ifttt_apikey" in config.keys():
        trigger_ifttt(event_name="snakemake_failed", apikey=config["ifttt_apikey"])
localrules: all, clean, copy_report_file, copy_report_file_without_resolution, report, supplemental_material
ruleorder: copy_report_file > copy_report_file_without_resolution

include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/uncertainty.smk"
include: "./rules/analyse.smk"
include: "./rules/impossible.smk"
include: "./rules/publication.smk"


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/logs/{resolution}/test-report.html".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/report.docx".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/report.pdf".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/supplemental.pdf".format(resolution=config["resolution"]["space"]),
        "build/output/national/uncertainty/time-diff.csv",
        "build/output/{resolution}/uncertainty/weather-diff-diff.txt".format(resolution=config["weather-uncertainty"]["resolution"]["space"])
        "build/output/{resolution}/figures/figure1.tiff".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/figures/figure2.tiff".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/figures/figure3.tiff".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/figures/figure4.tiff".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/figures/figure5.tiff".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/figures/figure6.tiff".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/figures/figure7.tiff".format(resolution=config["resolution"]["space"])


rule copy_report_file:
    message: "Copy file {input[0]} into dedicated report folder."
    input: "build/output/{resolution}/{filename}.{suffix}"
    wildcard_constraints: suffix = "((csv)|(png)|(svg)|(tiff))"
    output: "build/output/{resolution}/report/{filename}.{suffix}"
    shell: "cp {input} {output}"


rule copy_report_file_without_resolution:
    message: "Copy file {input[0]} into dedicated report folder."
    input: "build/output/{filename}.{suffix}"
    wildcard_constraints: suffix = "((csv)|(png)|(svg)|(tiff))"
    output: "build/output/{resolution}/report/{filename}.{suffix}"
    shell: "cp {input} {output}"


GENERAL_DOCUMENT_DEPENDENCIES = [
    "report/literature.bib",
    "report/report.css",
    "report/joule.csl",
    "report/template.html",
    "report/pandoc-metadata.yaml",
    "report/fonts/KlinicSlabBook.otf",
    "report/fonts/KlinicSlabBookIt.otf",
    "report/fonts/KlinicSlabMedium.otf",
    "report/fonts/KlinicSlabMediumIt.otf",
]


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--self-contained --css=report.css --template template.html --to html5"
    elif suffix == "pdf":
        return "--css=report.css --template template.html --pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        GENERAL_DOCUMENT_DEPENDENCIES,
        "report/report.md",
        "report/summaries.yaml",
        "build/output/{resolution}/report/total-sobol-diff.svg",
        "build/output/{resolution}/report/cost.svg",
        "build/output/{resolution}/report/map-cost.png",
        "build/output/{resolution}/report/map-energy.png",
        "build/output/{resolution}/report/composition.svg",
        "build/output/{resolution}/report/composition-uncertainty.png",
        "build/output/{resolution}/report/cost-uncertainty.png",
    output: "build/output/{resolution}/report.{suffix}"
    params: options = pandoc_options
    conda: "envs/pdf.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build/output/{wildcards.resolution}/report .
        {PANDOC} report.md {params.options} \
        --metadata-file=pandoc-metadata.yaml --metadata-file=summaries.yaml \
        -o ../build/output/{wildcards.resolution}/report.{wildcards.suffix}
        """


rule supplemental_material:
    message: "Compile the supplemental material."
    input:
        GENERAL_DOCUMENT_DEPENDENCIES,
        "report/supplemental.md",
        "report/supplemental.css",
        "report/biofuel-feedstocks.csv",
        "build/output/{resolution}/report/cost-special-cases.svg",
        "build/output/{resolution}/report/composition-all.svg",
        "build/output/{resolution}/report/total-sobol-all.svg",
        "build/output/{resolution}/report/first-sobol-all.svg",
        "build/output/{resolution}/report/total-minus-first-sobol-all.svg",
        "build/output/{resolution}/report/overview-scenario-results-1.csv",
        "build/output/{resolution}/report/overview-scenario-results-2.csv",
        "build/output/{resolution}/report/overview-uncertain-parameters.csv",
        "build/output/{resolution}/report/overview-cost-assumptions.csv",
        "build/output/{resolution}/report/network.png",
        "build/output/{resolution}/report/timeseries.svg",
        "build/output/{resolution}/report/generation-shares.svg",
        "build/output/{resolution}/report/bioenergy-use.svg"
    params: options = pandoc_options
    output: "build/output/{resolution}/supplemental.{suffix}"
    conda: "envs/pdf.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build/output/{wildcards.resolution}/report .
        {PANDOC} supplemental.md {params.options} --css=supplemental.css \
        --metadata-file=pandoc-metadata.yaml \
        -o ../build/output/{wildcards.resolution}/supplemental.{wildcards.suffix}
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """


def trigger_ifttt(event_name, apikey):
    import requests
    response = requests.post(
            f'https://maker.ifttt.com/trigger/{event_name}/with/key/{apikey}',
            data={"value1": "geographic-scale"}
    )
