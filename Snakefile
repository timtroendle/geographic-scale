PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

configfile: "./config/default.yaml"
wildcard_constraints:
    resolution = "((continental)|(national)|(regional))", # supported spatial resolutions
    scenario = "({})".format("|".join([f"({scenario})" for scenario in config["scenarios"]]))
onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'geographical-scale succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'geographical-scale crashed' {config[email]}")
localrules: all, clean, copy_report_file, report

include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/analyse.smk"
include: "./rules/impossible.smk"
include: "./rules/uncertainty.smk"


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/logs/{resolution}/test-report.html".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/report.html".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/aggregation.nc".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/uncertainty/weather-diff-diff.txt".format(resolution=config["weather-uncertainty"]["resolution"]["space"])


rule copy_report_file:
    message: "Copy file {input[0]} into dedicated report folder."
    input: "build/output/{resolution}/{filename}.{suffix}"
    wildcard_constraints: suffix = "((csv)|(png))"
    output: "build/output/{resolution}/report/{filename}.{suffix}"
    shell: "cp {input} {output}"


REPORT_DEPENDENCIES = [
    "report/report.md",
    "report/literature.bib",
    "report/biofuel-feedstocks.csv",
    "report/concept.md",
    "report/pandoc-metadata.yml",
    "report/report.css",
    "report/pnas.csl",
    "report/fonts/KlinicSlabBook.otf",
    "report/fonts/KlinicSlabBookIt.otf",
    "report/fonts/KlinicSlabMedium.otf",
    "report/fonts/KlinicSlabMediumIt.otf",
    "build/output/{resolution}/report/scenario-space.png",
    "build/output/{resolution}/report/map.png",
    "build/output/{resolution}/report/flows.png",
    "build/output/{resolution}/report/overview-scenario-results-1.csv",
    "build/output/{resolution}/report/overview-scenario-results-2.csv",
    "build/output/{resolution}/report/overview-cost-assumptions.csv"
]


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--self-contained --css=report.css --to html5"
    elif suffix == "pdf":
        return "--css=report.css --pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input: REPORT_DEPENDENCIES
    output: "build/output/{resolution}/report.{suffix}"
    params: options = pandoc_options
    conda: "envs/pdf.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build/output/{wildcards.resolution}/report .
        {PANDOC} report.md pandoc-metadata.yml {params.options} \
        -o ../build/output/{wildcards.resolution}/report.{wildcards.suffix}
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """
