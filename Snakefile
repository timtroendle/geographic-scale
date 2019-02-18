PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"


rule all:
    message: "Run entire analysis and compile report."
    input: "build/report.html"


rule scenario_space:
    message: "Plot sketch of scenario space."
    input: src="src/scenarios.py"
    output: "build/scenario-space.png"
    script: "src/scenarios.py"


rule map:
    message: "Create map of results."
    input:
        src = "src/map.py",
        continental_shape = "data/units-continental.geojson",
        national_shapes = "data/units-national.geojson",
        regional_shapes = "data/units-regional.geojson"
    output: "build/map.png"
    conda: "envs/geo.yaml"
    script: "src/map.py"


rule sensitivity_heatmap:
    message: "Plot sketch of sensitivity heatmap."
    input: src="src/sensitivity.py"
    output: "build/parameter-sensitivity.png"
    script: "src/sensitivity.py"


rule report:
    message: "Compile report."
    input:
        "report/literature.bib",
        "report/main.md",
        "report/pandoc-metadata.yml",
        "report/report.css",
        rules.scenario_space.output,
        rules.sensitivity_heatmap.output,
        rules.map.output
    output:
        "build/report.html"
    shell:
        """
        cd ./report
        {PANDOC} --self-contained --css=report.css main.md pandoc-metadata.yml \
        --to html5 -o ../build/report.html
        """


rule pdf_report:
    message: "Compile PDF report."
    input: rules.report.output # to avoid repeating all dependencies
    output: "build/report.pdf"
    conda: "envs/pdf.yaml"
    shell:
        """
        cd ./report
        {PANDOC} --css=report.css main.md pandoc-metadata.yml --pdf-engine weasyprint \
        -o ../build/report.pdf
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
