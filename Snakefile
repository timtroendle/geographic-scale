PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/analyse.smk"
include: "./rules/concept.smk"
localrules: all, clean

onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'ose-model-comparison succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'ose-model-comparison crashed' {config[email]}")


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/logs/test-report.html"


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """
