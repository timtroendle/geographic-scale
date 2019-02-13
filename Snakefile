PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"


rule all:
    message: "Run entire analysis and compile report."
    input: "build/report.html"


rule scenario_space:
    message: "Plot sketch of scenario space."
    input: src="src/scenarios.py"
    output: "build/scenario-space.png"
    script: "src/scenarios.py"


rule report:
    message: "Compile report."
    input:
        "report/literature.bib",
        "report/main.md",
        "report/pandoc-metadata.yml",
        rules.scenario_space.output
    output:
        "build/report.html"
    shell:
        """
        cd ./report
        {PANDOC} --self-contained --css=report.css main.md pandoc-metadata.yml \
        --to html5 -o ../build/report.html
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """


rule test:
    shell:
        "py.test"
