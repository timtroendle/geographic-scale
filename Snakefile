PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

configfile: "./config/default.yaml"
wildcard_constraints:
    resolution = "((continental)|(national)|(regional))", # supported spatial resolutions
    scenario = "({})".format("|".join([f"({scenario})" for scenario in config["scenarios"]])),
    plot_suffix = "((png)|(svg))"
onstart:
    shell("mkdir -p build/logs/uncertainty")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'geographical-scale succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'geographical-scale crashed' {config[email]}")
localrules: all, clean, copy_report_file, copy_report_file_without_resolution, report, supplementary_material
ruleorder: copy_report_file > copy_report_file_without_resolution

include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/uncertainty.smk"
include: "./rules/analyse.smk"
include: "./rules/impossible.smk"


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/logs/{resolution}/test-report.html".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/report.pdf".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/supplementary.pdf".format(resolution=config["resolution"]["space"]),
        "build/output/{resolution}/uncertainty/weather-diff-diff.txt".format(resolution=config["weather-uncertainty"]["resolution"]["space"])


rule copy_report_file:
    message: "Copy file {input[0]} into dedicated report folder."
    input: "build/output/{resolution}/{filename}.{suffix}"
    wildcard_constraints: suffix = "((csv)|(png)|(svg))"
    output: "build/output/{resolution}/report/{filename}.{suffix}"
    shell: "cp {input} {output}"


rule copy_report_file_without_resolution:
    message: "Copy file {input[0]} into dedicated report folder."
    input: "build/output/{filename}.{suffix}"
    wildcard_constraints: suffix = "((csv)|(png)|(svg))"
    output: "build/output/{resolution}/report/{filename}.{suffix}"
    shell: "cp {input} {output}"


GENERAL_DOCUMENT_DEPENDENCIES = [
    "report/literature.bib",
    "report/report.css",
    "report/joule.csl",
    "report/template.html",
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
        "report/pandoc-metadata.yml",
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
        {PANDOC} report.md pandoc-metadata.yml {params.options} \
        -o ../build/output/{wildcards.resolution}/report.{wildcards.suffix}
        """


rule supplementary_material:
    message: "Compile the supplementary material."
    input:
        GENERAL_DOCUMENT_DEPENDENCIES,
        "report/supplementary.md",
        "report/biofuel-feedstocks.csv",
        "build/output/{resolution}/report/cost-special-cases.svg",
        "build/output/{resolution}/report/total-sobol-all.svg",
        "build/output/{resolution}/report/first-sobol-all.svg",
        "build/output/{resolution}/report/total-minus-first-sobol-all.svg",
        "build/output/{resolution}/report/overview-scenario-results-1.csv",
        "build/output/{resolution}/report/overview-scenario-results-2.csv",
        "build/output/{resolution}/report/overview-uncertain-parameters.csv",
        "build/output/{resolution}/report/overview-cost-assumptions.csv",
        "build/output/{resolution}/report/network.png",
        "build/output/{resolution}/report/timeseries.svg"
    params: options = pandoc_options
    output: "build/output/{resolution}/supplementary.{suffix}"
    conda: "envs/pdf.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build/output/{wildcards.resolution}/report .
        {PANDOC} supplementary.md {params.options} --table-of-contents \
        -o ../build/output/{wildcards.resolution}/supplementary.{wildcards.suffix}
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """
