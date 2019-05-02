PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/analyse.smk"
include: "./rules/uncertainty.smk"
localrules: all, clean

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
        "build/output/{}/scenario-space.png".format(config["resolution"]),
        "build/output/{}/map.png".format(config["resolution"])


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """
