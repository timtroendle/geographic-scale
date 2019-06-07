PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/analyse.smk"
include: "./rules/uncertainty.smk"
localrules: all, clean, copy_report_file, report, pdf_report, docx_report

wildcard_constraints:
        resolution = "((national)|(regional))" # can run on national or regional spatial resolution

onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'geographical-scale succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'geographical-scale crashed' {config[email]}")


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/logs/test-report.html",
        "build/report.html"


rule copy_report_file:
    message: "Copy file {input[0]} into dedicated report folder."
    input: "build/output/{}/{{filename}}.{{suffix}}".format(config["resolution"]["space"])
    wildcard_constraints: suffix = "((csv)|(png))"
    output: "build/output/report/{filename}.{suffix}"
    shell: "cp {input} {output}"

REPORT_DEPENDENCIES = [
    "report/report.md",
    "report/literature.bib",
    "report/concept.md",
    "report/pandoc-metadata.yml",
    "report/report.css",
    "report/pnas.csl",
    "report/fonts/KlinicSlabBook.otf",
    "report/fonts/KlinicSlabBookIt.otf",
    "report/fonts/KlinicSlabMedium.otf",
    "report/fonts/KlinicSlabMediumIt.otf",
    "build/output/report/scenario-space.png",
    "build/output/report/map.png",
    "build/output/report/overview-scenario-results.csv",
    "build/output/report/overview-cost-assumptions.csv"
]

rule report:
    message: "Compile report."
    input: REPORT_DEPENDENCIES
    output:
        "build/report.html"
    conda: "envs/pdf.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} --self-contained --css=report.css report.md pandoc-metadata.yml \
        --to html5 -o ../build/report.html
        """


rule pdf_report:
    message: "Compile PDF report."
    input: REPORT_DEPENDENCIES
    output: "build/report.pdf"
    conda: "envs/pdf.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} --css=report.css report.md pandoc-metadata.yml --pdf-engine weasyprint \
        -o ../build/report.pdf
        """


rule docx_report:
    message: "Compile DOCX report."
    input: REPORT_DEPENDENCIES
    output: "build/report.docx"
    conda: "envs/pdf.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} report.md pandoc-metadata.yml -o ../build/report.docx
        """

rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """
